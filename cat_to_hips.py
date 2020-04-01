#!/usr/bin/env python3

# steps to convert catalogue to HiPS catalogue:

# Output:
#  properties: basic properties
#  metadata.xml: different columns (plus designated ra, dec)
#  Moc.fits: different regions in output
#  Moc.json: different encoding of above
#   Norder...

import os.path
import datetime
import json

import numpy as N

from astropy.io import fits
from astropy_healpix import HEALPix
from astropy import units as u

builder = 'cat_to_hips.py'

def make_dir(dirname):
    try:
        os.mkdir(dirname)
    except OSError:
        pass

def getDate():
    return datetime.datetime.now().strftime('%Y-%m-%dT%H:%M%Z')

class Column:
    def __init__(self, name, data, metaattrs={}):
        self.name = name
        self.data = data
        self.metaattrs = metaattrs

metahdr='''<?xml version='1.0'?>
<VOTABLE version="1.2"
 xmlns="http://www.ivoa.net/xml/VOTable/v1.2">
<RESOURCE>
<TABLE>
'''

metafoot = '''<DATA>
<TABLEDATA>
</TABLEDATA>
</DATA>
</TABLE>
</RESOURCE>
</VOTABLE>
'''

def writeMetadata(outfilename, columns):
    """Write metadata file.

    TODO: Use XML library
    """

    with open(outfilename, 'w') as fout:
        fout.write(metahdr)
        for col in columns:
            out = '<FIELD name="%s"' % col.name
            for attr in sorted(col.metaattrs):
                out = out + ' %s="%s"' % (attr, col.metaattrs[attr])
            out += '/>\n'
            fout.write(out)
        fout.write(metafoot)

def writeProperties(outdir, title, columns, norder):
    outnow = getDate()
    with open(os.path.join(outdir, 'properties'), 'w') as fout:
        name = outdir.rstrip('/').split('/')[-1]
        fout.write('publisher_did = ivo://PRIVATE_USER/%s\n' % name)
        fout.write('dataproduct_type = catalog\n')
        fout.write('obs_title = %s\n' % title)
        fout.write('hips_service_url = %s\n' % name)
        fout.write('hips_builder = %s\n' % builder)
        fout.write('hips_release_data = %s\n' % outnow)
        fout.write('hips_frame = equatorial\n')
        fout.write('hips_cat_nrows = %i\n' % len(columns[0].data))
        fout.write('hips_order = %i\n' % norder)
        fout.write('hips_tile_format = tsv\n')
        fout.write('hips_initial_ra = 0.0\n')
        fout.write('hips_initial_dec = +0.0\n')
        fout.write('hips_status = public master unclonable\n')
        fout.write('hips_version = 1.4\n')
        fout.write('# Deprecated but still in use\n')
        fout.write('label = %s\n' % title)
        fout.write('coordsys = C\n')

def makeMoc(outdir, ras, decs, maxorder=7):
    """Write Moc coverage files."""

    print('Constructing MOC')

    ra_deg = ras * u.deg
    dec_deg = decs * u.deg

    # work out coverage of pixels
    hp = HEALPix(2**maxorder, order='nested')
    pix = hp.lonlat_to_healpix(ra_deg, dec_deg)
    pix = list(N.unique(pix))

    # combine sets of four pixels to lower orders
    order = maxorder
    out = []
    while order>0:
        nextpix = []
        i = 0
        while i < len(pix)-3:
            vi = pix[i]
            if vi%4==0 and pix[i+1]==vi+1 and pix[i+2]==vi+2 and pix[i+3]==vi+3:
                nextpix.append(vi//4)
                pix = pix[:i]+pix[i+4:]
            else:
                i += 1
        out.insert(0, (order, N.array(pix)))
        pix = nextpix
        if len(pix)==0:
            break
        order -= 1

    if len(pix) > 0:
        out.insert(0, (order, N.array(pix)))

    # write json
    out_json = {}
    for order, pixels in out:
        out_json[str(order)] = [int(i) for i in pixels]

    with open(os.path.join(outdir, 'Moc.json'), 'w') as fout:
        fout.write('#MOCORDER %i\n' % maxorder)
        json.dump(out_json, fout)

    # now write in FITS format
    fits_pixels = []
    for order, pixels in out:
        fits_pixels.append(4 * 4**order + pixels)
    fits_pixels = N.concatenate(fits_pixels)
    fits_pixels.sort()
    print(' Included', len(fits_pixels), 'healpix pixels, up to level', order)

    cols = fits.ColDefs([
        fits.Column(name='NPIX', format='J', array=fits_pixels)])    
    hdu = fits.BinTableHDU.from_columns(cols)
    hdr = hdu.header
    hdr['PIXTYPE'] = 'HEALPIX'
    hdr['ORDERING'] = 'NUNIQ'
    hdr['COORDSYS'] = 'C'
    hdr['MOCORDER'] = maxorder
    hdr['MOCTOOL'] = builder
    hdr['DATE'] = getDate()

    hdulist = fits.HDUList([fits.PrimaryHDU(), hdu])
    hdulist.writeto(os.path.join(outdir, 'Moc.fits'), overwrite=True)

def split_catalog(ra_v, dec_v, score_v, aim=100, minorder=1):
    """Do splitting into orders.

    ra_v, dec_v: ra and dec
    score_v: score (higher values more likely to be shown first.
    """

    print('Splitting catalog into orders')

    tot_num = len(ra_v)
    splits = N.zeros(tot_num, dtype=N.int32)-1

    # lower threshold by factor ratio until the median number of
    # objects in cels is greater than aim
    ratio = 0.99
    score = N.max(score_v)*ratio
    tot_sel = 0

    order = minorder
    while True:
        hp = HEALPix(2**order, order='nested')
        sel = (splits<0) & (score_v>=score)
        num_sel = N.sum(sel)
        if num_sel > 0:
            pix = hp.lonlat_to_healpix(ra_v[sel]*u.deg, dec_v[sel]*u.deg)
            pix, cts = N.unique(pix, return_counts=True)
            medcts = N.median(cts[cts>0])

        if (num_sel==0 or len(pix) < 4 or medcts < aim) and tot_sel+num_sel != tot_num:
            score *= ratio
        else:

            print(N.percentile(cts[cts>0], [5, 50, 95]))

            splits[sel] = order
            tot_sel += num_sel
            print(' setting order=%i using score=%g (%i objects)' % (
                order, score, num_sel))
            if N.sum(splits<0) == 0:
                break
            order += 1

    return order, splits

def write_cat(outfilename, hdr, columns, sel_idx):
    with open(outfilename, 'w') as fout:
        fout.write(hdr)
        for idx in sel_idx:
            line = '\t'.join([str(c.data[idx]) for c in columns]) + '\n'
            fout.write(line)

def write_order(outdir, ra_v, dec_v, columns, order, order_splits):
    print('Writing order', order)
    rootorder = os.path.join(outdir, 'Norder%i' % order)
    make_dir(rootorder)

    colnames = [c.name for c in columns]
    hdr = '\t'.join(colnames) + '\n'

    # indices of which sources to include
    ordersel = N.where(order_splits == order)[0]

    # write all sky file
    if order<=3:
        allskyfile = os.path.join(rootorder, 'Allsky.tsv')
        write_cat(allskyfile, hdr, columns, ordersel)

    # split into pixels
    hp = HEALPix(2**order, order='nested')
    pixels = hp.lonlat_to_healpix(ra_v[ordersel]*u.deg, dec_v[ordersel]*u.deg)

    madedir = set()

    # iterate over each pixel for the order
    for pix in range(0, 12*4**order):
        diridx = (pix // 10000)*10000
        thedir = os.path.join(rootorder, 'Dir%i' % diridx)
        if diridx not in madedir:
            print(' Making', thedir)
            make_dir(thedir)
            madedir.add(diridx)
        outpixfile = os.path.join(thedir, 'Npix%i.tsv' % pix)

        # lookup subset in subset
        idxs = ordersel[N.where(pixels==pix)[0]]
        write_cat(outpixfile, hdr, columns, idxs)

def catalog_to_hips(outdir, columns,
                    title='My Cat',
                    ra_col='RA', dec_col='DEC', score_col='SCORE'):

    """Convert catalog into HiPS format.

    outdir: output directory
    columns: list of Column objects
    title: catalog title
    ra_col: name of column for RA (deg)
    dec_col: name of column for Dec (deg)
    score_col: name of column for (increasing) score
    """

    print('Converting catalog with %i sources' % len(columns[0].data))

    make_dir(outdir)

    colmap = {col.name: col for col in columns}

    # specify position columns for metadata
    colmap[ra_col].metaattrs['ucd'] = 'pos.eq.ra;meta.main'
    colmap[ra_col].metaattrs['unit'] = 'degree'
    colmap[ra_col].metaattrs['datatype'] = 'double'
    colmap[dec_col].metaattrs['ucd'] = 'pos.eq.dec;meta.main'
    colmap[dec_col].metaattrs['unit'] = 'degree'
    colmap[dec_col].metaattrs['datatype'] = 'double'

    writeMetadata(os.path.join(outdir, 'Metadata.xml'), columns)
    writeMetadata(os.path.join(outdir, 'metadata.xml'), columns)

    ra_v = colmap[ra_col].data
    dec_v = colmap[dec_col].data
    score_v = colmap[score_col].data

    norder, order_splits = split_catalog(ra_v, dec_v, score_v, minorder=1)

    for order in range(N.min(order_splits), max(3, norder)+1):
        write_order(outdir, ra_v, dec_v, columns, order, order_splits)

    makeMoc(outdir, ra_v, dec_v, maxorder=7)
    writeProperties(outdir, title, columns, norder=norder)

def main():
    f = fits.open('cat2rxs.fits', 'readonly')

    ra_deg = f[1].data['RA_DEG']
    dec_deg = f[1].data['DEC_DEG']
    exi_ml = f[1].data['EXI_ML']

    # sel = N.arange(len(ra_deg))
    # N.random.shuffle(sel)
    # sel = sel[:1000]
    # ra_deg = ra_deg[sel]
    # dec_deg = dec_deg[sel]
    # exi_ml = exi_ml[sel]

    cols = [
        Column('RA_DEG', ra_deg),
        Column('DEC_DEG', dec_deg),
        Column('EXI_ML', exi_ml, metaattrs={'datatype': 'float'}),
        ]

    catalog_to_hips('mycat', cols, ra_col='RA_DEG', dec_col='DEC_DEG', score_col='EXI_ML')

if __name__ == '__main__':
    main()

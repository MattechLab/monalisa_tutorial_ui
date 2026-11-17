"""Read-only verification of the six selected 13-November scanner exports."""
from pathlib import Path
import gzip
import hashlib
import json
import struct
import sys
import numpy as np

root=Path(sys.argv[1]).resolve()
expected={702:('00016','15689',6.7,2.875,12),
          706:('00017','15690',13.,2.875,12),
          707:('00018','15691',13.,4.,12),
          705:('00019','15692',6.7,2.875,20),
          708:('00020','15693',13.,2.875,20),
          709:('00021','15694',13.,4.,20)}
actual={p.name for p in root.iterdir() if p.is_dir() and p.name.startswith('seq')}
assert actual=={f'seq{sid}' for sid in expected},actual
rows=[]
for sid,(mid,fid,tr,te,fa) in expected.items():
    folder=root/f'seq{sid}'
    report=json.loads((folder/'reconstruction.json').read_text())
    nii=folder/'xrms_normalized.nii.gz'; metadata=report['metadata']
    assert hashlib.md5(nii.read_bytes()).hexdigest()==report['niftiMD5']
    with gzip.open(nii,'rb') as stream:
        header=stream.read(352)
        data=np.frombuffer(stream.read(),dtype='<f4')
    assert struct.unpack_from('<i',header,0)[0]==348
    assert struct.unpack_from('<4h',header,40)==(3,240,240,240)
    assert struct.unpack_from('<h',header,70)[0]==16
    assert struct.unpack_from('<3f',header,80)==(1.,1.,1.)
    assert struct.unpack_from('<f',header,108)[0]==352.
    assert struct.unpack_from('<2h',header,252)==(0,2)
    for offset,values in [(280,(1.,0.,0.,-120.)),(296,(0.,1.,0.,-120.)),(312,(0.,0.,1.,-120.))]:
        assert struct.unpack_from('<4f',header,offset)==values
    assert data.size==240**3 and np.isfinite(data).all()
    assert data.min()>=0 and data.max()==1
    assert report['matrixSize']==[240]*3 and report['spacingMm']==[1]*3
    assert report['normalizedArrayExactReadback'] and report['relativeRawIntensityRestorationError']<1e-7
    assert report['normalizationDivisor']>0
    assert metadata['nShotOff']==50 and metadata['physicalCoils']==32
    assert metadata['N']==480 and metadata['nSeg']==44 and metadata['nShot']==500
    assert metadata['retainedShots']==450 and metadata['imagingReadouts']==19350
    assert metadata['kernel']=='gauss' and metadata['kernelWidth']==3
    assert metadata['kernelParameters']==[.61,10]
    assert not metadata['signalConjugated'] and not metadata['additionalADCPhaseCorrection']
    assert not metadata['gradientDelayCorrection'] and not metadata['coilCompression']
    raw=Path(metadata['rawFile']); sequence=Path(metadata['seqFile'])
    assert raw.name==f'meas_MID{mid}_FID{fid}_seq{sid}.dat'
    assert raw.parent==Path('/Volumes/Extreme SSD/261113')
    assert sequence.parent==raw.parent and f'seq{sid}_' in sequence.name
    assert raw.stat().st_size==metadata['rawStat']['bytes']
    assert hashlib.md5(sequence.read_bytes()).hexdigest()==metadata['sequenceMD5']
    h=metadata['sequenceParams']['header']
    assert abs(h['TR_s']*1000-tr)<1e-8 and abs(h['TE_s']*1000-te)<1e-8
    assert metadata['sequenceFileFlipAngleDeg']==h['rf_flip_deg']==fa/2
    assert metadata['reportedFlipAngleDeg']==fa
    assert not (folder/'reconstruction.mat').exists()
    assert not (folder/'xrms_normalized.nii').exists()
    for name in ['reconstruction.log','xrms_preview.png']:
        assert (folder/name).stat().st_size>0
    rows.append({'sequence':sid,'nifti':str(nii),'TR_ms':tr,'TE_ms':te,
                 'reported_FA_deg':fa,'sequence_file_FA_deg':fa/2,
                 'MiB':nii.stat().st_size/1048576,
                 'normalizationDivisor':report['normalizationDivisor'],
                 'relative_restoration_error':report['relativeRawIntensityRestorationError'],
                 'verified':True})
    print(f'seq{sid}: VERIFIED; TR/TE {tr:g}/{te:g} ms; reported FA {fa}°; {rows[-1]["MiB"]:.2f} MiB')
result={'batch':str(root),'checks':'Six exact source pairings; NIfTI dimensions/float32/spacing/affine/payload/MD5/gzip; finite max-normalized values; acquisition/reconstruction metadata; corrected reporting angles; unchanged sequence hashes and raw sizes; compact file inventory. MATLAB exporter exact read-back passed.',
        'total_MiB':sum(r['MiB'] for r in rows),'results':rows}
(Path(__file__).resolve().parent/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(f'All six verified; {result["total_MiB"]:.2f} MiB total.')

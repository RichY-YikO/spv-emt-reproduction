from pathlib import Path
import argparse
import hashlib
import json
import shutil
import sys
import tempfile

ROOT=Path(__file__).resolve().parent
CODE=ROOT/'code'
sys.path.insert(0,str(CODE))
FIGURES=[f'Fig{i}' for i in range(1,6)]+[f'FigS{i}' for i in range(1,10)]

def main():
    parser=argparse.ArgumentParser(description='Redraw final publication figures from the bundled caches without recomputing scientific results.')
    parser.add_argument('--figures',nargs='+',choices=FIGURES+['Fig1a'],default=FIGURES)
    parser.add_argument('--output-dir',type=Path,default=ROOT/'outputs')
    args=parser.parse_args()
    output=args.output_dir.resolve()
    for sub in ['pdf','png','validation/layout']:(output/sub).mkdir(parents=True,exist_ok=True)
    import reproduce_fig1_fig2
    from replot_publication import replot
    selected=list(dict.fromkeys(args.figures))
    with tempfile.TemporaryDirectory(prefix='spv-emt-redraw-') as folder:
        work=Path(folder)
        ordinary=[name for name in selected if name!='Fig1a']
        records=[]
        if ordinary:
            replot(work,ordinary)
            manifest=json.loads((work/'replot_manifest.json').read_text(encoding='utf-8'))
            for rel in manifest['last_redraw']:
                base=work/rel
                records.append((base.name,base))
        if 'Fig1a' in selected:
            import numpy as np
            from reproduce_fig1a_no_orientation import require_pearson_cache,draw
            cache=CODE/'fig1_fig2_results'/'S_final_pearson.npz'
            require_pearson_cache(cache)
            with np.load(cache) as z:draw(work,z['signed_maxabs_r'],z['extreme_null_lo'],z['extreme_null_hi'])
            records.append(('Fig1a',work/'Fig1a_no_orientation'))
        for name,base in records:
            for ext,kind in [('.pdf','pdf'),('.png','png'),('.layout.json','validation/layout')]:
                shutil.copy2(base.with_suffix(ext),output/kind/(name+ext))
    hashes={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in CODE.rglob('*.npz')}
    (output/'redraw_manifest.json').write_text(json.dumps({'selection':selected,'scientific_inference_recomputed':False,'cache_sha256':hashes},indent=2),encoding='utf-8')
    print('Final figures saved to '+str(output),flush=True)

if __name__=='__main__':
    main()

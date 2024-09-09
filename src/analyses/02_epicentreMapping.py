import numpy as np
import pickle
from enigmatoolbox.datasets import load_sc, load_fc


# Help function
def load_result():
    
    
# Main analysis
def main():

    print('---------------------------')
    print('Polygenic epicenter mapping')
    print('---------------------------')

    fc_ctx, _, fc_sctx, _ = load_fc()
    sc_ctx, _, sc_sctx, _ = load_sc()
    prs_map =

    print()
    print('Functional networks')
    print('---------------------------')
    # Cortical epicenters
    fc_ctx_r, fc_ctx_p = epicenter_mapping(prs_map, fc_ctx)
    # Subcortical epicenters
    fc_sctx_r, fc_sctx_p = epicenter_mapping(prs_map, fc_sctx)

    print()
    print('Structural networks')
    print('---------------------------')
    # Cortical epicenters
    sc_ctx_r, sc_ctx_p = epicenter_mapping(prs_map, sc_ctx)
    # Subcortical epicenters
    sc_sctx_r, sc_sctx_p = epicenter_mapping(prs_map, sc_sctx)

    print()
    print('Save results')
    print('---------------------------')
    np.savez('../../data/results/02_epicentreMapping/epicentre.npz',
             fc_ctx_r=fc_ctx_r,  fc_ctx_p=fc_ctx_p,
             fc_sctx_r=fc_sctx_r, fc_sctx_p=fc_sctx_p,
             sc_ctx_r=sc_ctx_r, sc_ctx_p=sc_ctx_p,
             sc_sctx_r=sc_sctx_r, sc_sctx_p=sc_sctx_p)


if __name__ == "__main__":
    main()

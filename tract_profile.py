"""
==========================
Plotting tract profiles
==========================

An example of tracking and segmenting two tracts, and plotting their tract
profiles for FA (calculated with DTI).

"""
import os.path as op
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import nibabel as nib
import dipy.data as dpd
from dipy.data import fetcher
from dipy.io.gradients import read_bvals_bvecs
from dipy.core.gradients import gradient_table

import AFQ.utils.streamlines as aus
import AFQ.data as afd
import AFQ.tractography as aft
import AFQ.registration as reg
import AFQ.dti as dti
import AFQ.segmentation as seg

plt.switch_backend('agg')

data_path = '/N/dc2/projects/lifebid/HCP7/108323/'
data_file = data_path + 'original_hcp_data/Diffusion_7T/' + 'data.nii.gz'
data_bvec = data_path + 'original_hcp_data/Diffusion_7T/' + 'bvecs'
data_bval = data_path + 'original_hcp_data/Diffusion_7T/' + 'bvals'

#data_brainmask = data_path + 'original_hcp_data/Diffusion_7T/' + 'nodif_brain_mask.nii.gz'
#fs_path = '/N/dc2/projects/lifebid/HCP/Brent/7t_rerun/freesurfer/7t_108323'
#data_fs_seg = fs_path + '/mri/aparc+aseg.nii.gz'

"""
dpd.fetch_stanford_hardi()

hardi_dir = op.join(fetcher.dipy_home, "stanford_hardi")
hardi_fdata = op.join(hardi_dir, "HARDI150.nii.gz")
hardi_fbval = op.join(hardi_dir, "HARDI150.bval")
hardi_fbvec = op.join(hardi_dir, "HARDI150.bvec")
"""

img = nib.load(data_file)

print("Calculating DTI...")
if not op.exists('./dti_FA.nii.gz'):
    dti_params = dti.fit_dti(hardi_fdata, hardi_fbval, hardi_fbvec,
                             out_dir='.')
else:
    dti_params = {'FA': './dti_FA.nii.gz',
                  'params': './dti_params.nii.gz'}

tg = nib.streamlines.load('./dti_streamlines.tck').tractogram
streamlines = tg.apply_affine(np.linalg.inv(img.affine)).streamlines

# Use only a small portion of the streamlines, for expedience:
streamlines = streamlines[::100]

templates = afd.read_templates()
bundle_names = ["CST", "ILF"]

bundles = {}
for name in bundle_names:
    for hemi in ['_R', '_L']:
        bundles[name + hemi] = {'ROIs': [templates[name + '_roi1' + hemi],
                                         templates[name + '_roi1' + hemi]],
                                'rules': [True, True]}


print("Registering to template...")
MNI_T2_img = dpd.read_mni_template()
bvals, bvecs = read_bvals_bvecs()
gtab = gradient_table(bvals, bvecs, b0_threshold=100)
mapping = reg.syn_register_dwi(data_file, gtab)
reg.write_mapping(mapping, './mapping.nii.gz')


print("Segmenting fiber groups...")
fiber_groups = seg.segment(hardi_fdata,
                           hardi_fbval,
                           hardi_fbvec,
                           streamlines,
                           bundles,
                           reg_template=MNI_T2_img,
                           mapping=mapping,
                           as_generator=False,
                           affine=img.affine)


FA_img = nib.load(dti_params['FA'])
FA_data = FA_img.get_data()

print("Extracting tract profiles...")
for bundle in bundles:
    fig, ax = plt.subplots(1)
    profile = seg.calculate_tract_profile(FA_data, fiber_groups[bundle])
    ax.plot(profile)
    ax.set_title(bundle)
plt.savefig('bundle.png')
plt.show()
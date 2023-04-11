# -*- coding: utf-8 -*-
import pandas as pd
import copy
import nibabel as nib
import os
import argparse
import json
from multiprocessing import Pool as Pool
from nipype.interfaces import afni
import nipype.interfaces.spm as spm

def cli_parser():
	parser = argparse.ArgumentParser(description = __doc__)
	parser.add_argument('-toolbox_dir', help = 'The toolbox directory.')
	parser.add_argument('-input', help = 'The abs dir of input directory.')
	parser.add_argument('-output', help = 'The abs dir of output directory.')
#	parser.add_argument('-sublist', help = 'The path to the subjects list file.')
	parser.add_argument('-num', help = 'The number of which subject you want begin with.')
	parser.add_argument('-n', help = 'The total number of subjects.')
	parser.add_argument('-res', default = 3,help = 'The resolution of the template for function image normalization. Default 3 mm')
	parser.add_argument('-func_id', help = 'Select a specific task to be processed')
	parser.add_argument('-SliceTiming', default = 'NoSliceTiming', help = 'The path to Slice Order file or NoSliceTiming. Default NoSliceTiming.')
	parser.add_argument('-rm_TP', help = 'The number of the first time point you want to remove.')
	parser.add_argument('-fieldmap', default = 0, help = 'Field Distortions Correction, 1:True; 0: False. Default 1.')
	parser.add_argument('-mc_method', default = 'fsl', help = 'The method of Motion Correction, "fsl", "afni","spm","None".Default fsl.')
	parser.add_argument('-fwhm',default = 6, help = 'Width of gaussian kernel (mm).Default 6 mm.')
	parser.add_argument('-denoise',nargs = 4, default= ['ICA',0,0,0], help = 'The method of denoise.["no, regress, or ICA", motion: 6,12,24, \
		CSF: 0:mean; -1: mean,t-1 and quadratic form; positive num: PC num, WM: Same as CSF. Default ICA denoise')
	parser.add_argument('-despike', default = 0,help = 'Wavelet Despike, 1:True; 0: False.Default 0.')
	parser.add_argument('-filter',nargs = 2, default= [0,0], help = 'Filter. [highpass,lowpass], if you want only highpass filter, set lowpass = 0; \
		if you jsut want to detrend, set both lowpass and highpass = 0')

	return parser

parser = cli_parser()
args = parser.parse_args()

#Parameters Statement
toolbox_dir = os.path.join(args.toolbox_dir,'')
input_directory = os.path.join(args.input,'')
output_directory = os.path.join(args.output,'')
sublist = os.listdir(args.input)
sublist.sort()
#~ f = open(args.sublist)
#~ sublist = f.read()
#~ sublist = sublist.split('\n')
num = int(args.num)
n = int(args.n)
Res = str(args.res)
func_id = str(args.func_id)
SliceOrder = str(args.SliceTiming)
first_TP = int(args.rm_TP)
fieldmap = bool(int(args.fieldmap))
mc_method = str(args.mc_method)
fwhm = float(args.fwhm)
denoise = list(args.denoise)
despike = bool(int(args.despike))
denoise_method = str(denoise[0])
motion = int(denoise[1])
CSF = int(denoise[2])
WM = int(denoise[3])
fil = list(args.filter)
highpass = float(fil[0])
lowpass = float(fil[1])

Standard_image_T1w = os.path.join(toolbox_dir,'MNI152/MNI152_T1_1mm_brain.nii.gz')
Standard_t1w = os.path.join(toolbox_dir,'MNI152/MNI152_T1_1mm.nii.gz')
Standard_mask = os.path.join(toolbox_dir,'MNI152/MNI152_T1_1mm_brain_mask.nii.gz')
Standard_image_func = os.path.join(toolbox_dir,'MNI152/MNI152_T1_'+Res+'mm_brain.nii.gz')
Standard_func_mask = os.path.join(toolbox_dir,'MNI152/MNI152_T1_'+Res+'mm_brain_mask.nii.gz')

os.system('mkdir '+output_directory)

def T1w_prep(input_dir,output_dir,sub):
	T1_directory = os.path.join(input_dir,'anat','')
	T1_files = os.listdir(T1_directory)
	for i in range(len(T1_files)):
		if T1_files[i].find('T1w.nii') != -1 :
			T1_name = T1_files[i]
			T1_file = os.path.join(T1_directory,T1_files[i])

	os.system('mkdir '+output_dir+'tmp/T1')
	os.system('mkdir '+output_dir+'T1')

    #preprocess T1 image
    #Crop
	if os.path.exists(output_dir+'tmp/T1/crop_'+T1_name) == 0:
		os.system('robustfov -i '+T1_file+' -r '+output_dir+'tmp/T1/crop_'+T1_name)

    #Skullstrip
	if os.path.exists(output_dir+'tmp/T1/BEBrainExtractionBrain.nii.gz') == 0:
		os.system('antsBrainExtraction.sh -d 3 -a '+output_dir+'tmp/T1/crop_'+T1_name+' -e '+Standard_t1w+' -m '+Standard_mask+' -o '+output_dir+'tmp/T1/BE')

    #Registration
	if os.path.exists(output_dir+'T1/T1w2std_'+T1_name) == 0:
		os.system('antsRegistrationSyN.sh -d 3 -f '+Standard_image_T1w+' -m '+output_dir+'tmp/T1/BEBrainExtractionBrain.nii.gz  -o ' \
			+output_dir+'tmp/T1/rega2t -n 20')
		os.system('antsApplyTransforms -d 3 -i '+output_dir+'tmp/T1/crop_'+T1_name+' -o ['+output_dir+'T1/T1w2std_Warp.nii.gz,1] -r '+Standard_t1w+ \
				' -t '+output_dir+'tmp/T1/rega2t1Warp.nii.gz -t '+output_dir+'tmp/T1/rega2t0GenericAffine.mat')
		os.system('antsApplyTransforms -d 3 -i '+output_dir+'tmp/T1/crop_'+T1_name+' -o '+output_dir+'T1/T1w2std_'+T1_name+' -r '+Standard_t1w+ \
				' -t '+output_dir+'T1/T1w2std_Warp.nii.gz')

    #Segmentation
	os.system('mkdir '+output_dir+'Seg')

	if os.path.exists(output_dir+'Seg/'+sub+'_CSF_mask.nii.gz') == 0 or os.path.exists(output_dir+'Seg/'+sub+'_WM_mask.nii.gz') == 0:
		os.system('fast -S 1 -t 1 -o '+output_dir+'Seg/seg -g -n 3 -b -I 10 '+output_dir+'tmp/T1/rega2tWarped.nii.gz')
		os.system('mv '+output_dir+'Seg/seg_seg_0.nii.gz '+output_dir+'Seg/'+sub+'_CSF_mask.nii.gz')
		os.system('mv '+output_dir+'Seg/seg_seg_2.nii.gz '+output_dir+'Seg/'+sub+'_WM_mask.nii.gz')
    
    # #Reslice the mask
	os.system('3dresample -dxyx '+Res+' '+Res+' '+Res+' -prefix '+output_dir+'Seg/'+sub+'_CSF_mask_'+Res+'mm.nii.gz -input '+output_dir+'Seg/'+sub+'_CSF_mask.nii.gz')
	os.system('3dresample -dxyx '+Res+' '+Res+' '+Res+' -prefix '+output_dir+'Seg/'+sub+'_WM_mask_'+Res+'mm.nii.gz -input '+output_dir+'Seg/'+sub+'_WM_mask.nii.gz')

def Prepare_fieldmap(input_dir,output_dir):
	fmap_directory = os.path.join(input_dir,'fmap','')
	fmap_files = os.listdir(fmap_directory)
	magnitude_file = list()
	phasediff_file = list()
	AP_file = list()
	PA_file = list()
	fieldmap_method = ''
	for i in range(len(fmap_files)):
		if fmap_files[i].find('magnitude1.nii') != -1 or fmap_files[i].find('magnitude.nii') != -1 :
			magnitude_name = fmap_files[i]
			magnitude_file = os.path.join(fmap_directory,fmap_files[i])

		if fmap_files[i].find('phasediff.nii') != -1 :
			phasediff_name = fmap_files[i]
			phasediff_file = os.path.join(fmap_directory,fmap_files[i])

		if fmap_files[i].find('epi.nii') != -1 and fmap_files[i].find('dir-AP') != -1 :
			AP_name = fmap_files[i]
			AP_file = os.path.join(fmap_directory,fmap_files[i])

		if fmap_files[i].find('epi.nii') != -1 and fmap_files[i].find('dir-PA') != -1 :
			PA_name = fmap_files[i]
			PA_file = os.path.join(fmap_directory,fmap_files[i])

	os.system('mkdir '+output_dir+'tmp/fmap')

	if magnitude_file and phasediff_file:
		mjson = os.path.join(fmap_directory,magnitude_name.replace('nii.gz','json'))
		pjson = os.path.join(fmap_directory,phasediff_name.replace('nii.gz','json'))
		with open(mjson,'r') as m:
			load_m = json.load(m)
		TE_m = load_m['EchoTime']
		with open(pjson,'r') as p:
			load_p = json.load(p)
		TE_p = load_p['EchoTime']
		dt = round(1000*(TE_p - TE_m),2)

		magnitude_brain = magnitude_name.replace('.nii.gz','_brain.nii.gz')
		magnitude_ero = magnitude_brain.replace('.nii.gz','_ero.nii.gz')
		os.system('bet '+magnitude_file+' '+output_dir+'tmp/fmap/'+magnitude_brain+' -f 0.4 -n -R')
		os.system('fslmaths '+output_dir+'tmp/fmap/'+magnitude_brain+' -ero '+output_dir+'tmp/fmap/'+magnitude_ero)
		os.system('fsl_prepare_fieldmap SIEMENS '+phasediff_file+' '+output_dir+'tmp/fmap/'+magnitude_ero+' '+output_dir+'tmp/fmap/fieldmap.nii.gz '+str(dt))

		fieldmap_method = 'Fugue'

	if AP_file and PA_file:
		APjson = os.path.join(fmap_directory,AP_name.replace('nii.gz','json'))
		PAjson = os.path.join(fmap_directory,PA_name.replace('nii.gz','json'))
		with open(APjson,'r') as AP:
			load_AP = json.load(AP)
		AP_readout = load_AP['TotalReadoutTime']
		with open(PAjson,'r') as PA:
			load_PA = json.load(PA)
		PA_readout = load_PA['TotalReadoutTime']

		os.system('printf "0 -1 0 '+str(AP_readout)+str('\n')+'0 1 0 '+str(PA_readout)+'" > '+output_dir+'tmp/fmap/acqparams.txt')
		os.system('fslmerge -t '+output_dir+'tmp/fmap/fieldmap.nii.gz '+AP_file+' '+PA_file)
		os.system('topup --imain='+output_dir+'tmp/fmap/fieldmap.nii.gz --datain='+output_dir+'tmp/fmap/acqparams.txt --out='+output_dir+'tmp/fmap/fieldmap_topup --iout='+output_dir+'tmp/fmap/fieldmap_unwarp')

		fieldmap_method = 'Topup'

	return fieldmap_method

def Motion_corr(fMRI,output_dir,name,method):
	if method == 'fsl':
		os.system('mcflirt -in '+fMRI+' -out '+output_dir+'tmp/func/mc_'+name+' -refvol 0 -plots')
		os.system('mv '+output_dir+'tmp/func/mc_'+name+'.par '+output_dir+'func/mc_'+name+'.txt')

	if method == 'afni':
		os.system('3dvolreg -Fourier -twopass -prefix '+output_dir+'tmp/func/mc_'+name+' -1Dfile '+output_dir+'func/mc_'+name+'.txt '+fMRI)

	if method == 'spm':
		os.system('gunzip '+fMRI)
		fMRI = fMRI.replace('.gz','')
		name = name.replace('.gz','')
		realign = spm.Realign()
		realign.inputs.in_files = fMRI
		realign.inputs.register_to_mean = False
		realign.run()
		os.system('mv '+output_dir+'tmp/func/rrm_'+name+' '+output_dir+'tmp/func/mc_'+name)
		os.system('gzip '+output_dir+'tmp/func/mc_'+name)
		os.system('mv '+output_dir+'tmp/func/rp* '+output_dir+'func/mc_'+name+'.txt')

def Fieldmap_corr(fMRI,output_dir,data_dir,name,method):
	if method == 'Fugue':
		fjson = os.path.join(data_dir,'func',name.replace('nii.gz','json'))
		with open(fjson,'r') as f:
			load_f = json.load(f)
		EchoSpace = round(load_f['EffectiveEchoSpacing'],6)
		os.system('fugue -i '+fMRI+' --dwell='+str(EchoSpace)+' --unwarpdir=y- --loadfmap='+output_dir+'tmp/fmap/fieldmap.nii.gz -u '+output_dir+'tmp/func/fieldcorr_'+name)

	if method == 'Topup':
		os.system('applytopup --imain='+fMRI+' --datain='+output_dir+'tmp/fmap/acqparams.txt --topup='+output_dir+'tmp/fmap/fieldmap_topup --inindex=1 --method=jac --out='+output_dir+'tmp/func/fieldcorr_'+name)

	out_name = 'fieldcorr_'+name

	return out_name

def Normalization(fMRI,output_dir,name,ref,Vol,TR):
	if os.path.exists(output_dir+'tmp/func/mean_'+name) == 0:
		os.system('antsMotionCorr -d 3 -a '+fMRI+' -o '+output_dir+'tmp/func/mean_'+name)
    #regirter function image to T1 image
	if os.path.exists(output_dir+'tmp/func/regf2aWarped.nii.gz') == 0:
		os.system('antsRegistrationSyN.sh -d 3 -f '+output_dir+'tmp/T1/BEBrainExtractionBrain.nii.gz -m '+output_dir+'tmp/func/mean_'+name+' -t "a" -o '+output_dir+'tmp/func/regf2a -n 20')

    #estimating the transformation
	if os.path.exists(output_dir+'tmp/func/func2standardWarp.nii.gz') == 0:
		os.system('antsApplyTransforms -d 3 -i '+output_dir+'tmp/func/mean_'+name+' -o '+output_dir+'func/example_f2std_'+name+' -r '+ref+' -t '+output_dir+'T1/T1w2std_Warp.nii.gz -t ' \
			+output_dir+'tmp/func/regf2a0GenericAffine.mat')
		os.system('antsApplyTransforms -d 3 -i '+output_dir+'tmp/func/mean_'+name+' -o ['+output_dir+'func/f2std_Warp_'+name+',1] -r '+ref+' -t '+output_dir+'T1/T1w2std_Warp.nii.gz -t ' \
			+output_dir+'tmp/func/regf2a0GenericAffine.mat')

	if os.path.exists(output_dir+'tmp/func/template_replicated.nii.gz') == 0:
		os.system('ImageMath 3 '+output_dir+'tmp/func/template_replicated.nii.gz ReplicateImage '+ref+' '+str(Vol)+' '+str(TR)+' 0')

	if os.path.exists(output_dir+'tmp/func/4Dfunc2standardWarp.nii.gz') == 0:
		os.system('ImageMath 3 '+output_dir+'tmp/func/4Dfunc2standardWarp.nii.gz ReplicateDisplacement '+output_dir+'func/f2std_Warp_'+name+' '+str(Vol)+' '+str(TR)+' 0')

	if os.path.exists(output_dir+'tmp/func/f2std_'+name) == 0:
		os.system('antsApplyTransforms -d 4 -o '+output_dir+'tmp/func/f2std_'+name+' -t '+output_dir+'tmp/func/4Dfunc2standardWarp.nii.gz -r '+output_dir+'tmp/func/template_replicated.nii.gz -i '+fMRI)

	if os.path.exists(output_dir+'tmp/func/f2std_'+name) == 1:
		os.system('rm '+output_dir+'tmp/func/regf2aWarped.nii.gz')
		os.system('rm '+output_dir+'tmp/func/template_replicated.nii.gz')
		os.system('rm '+output_dir+'tmp/func/4Dfunc2standardWarp.nii.gz')

def Smooth(fMRI,output_dir,name,fwhm):
	smooth = afni.preprocess.BlurToFWHM()
	smooth.inputs.in_file = fMRI
	smooth.inputs.fwhm = fwhm
	smooth.inputs.outputtype = 'NIFTI_GZ'
	smooth.inputs.out_file = output_dir+'tmp/func/smooth_'+name
	smooth.inputs.num_threads = 20
	smooth.run()

	smooth_name = 'smooth_'+name

	return smooth_name

def prep_main(sub):
	data_directory = os.path.join(input_directory,sub)
	output = os.path.join(output_directory,sub,'')
	# if os.path.exists(output+'tmp') != 0:
	# 	os.system('rm -r '+output)

	os.system('mkdir '+output)
	os.system('mkdir '+output+'tmp')

	T1w_prep(data_directory,output,sub)

	if fieldmap:
		fieldmap_method = Prepare_fieldmap(data_directory,output)

	func_directory = os.path.join(data_directory,'func','')
	func_files = os.listdir(func_directory)
	os.system('mkdir '+output+'func')
	
	index = []
	for j in range(len(func_files)):
		if func_files[j].find('.nii') != -1 and func_files[j].find(func_id) != -1:
			index.append(func_files[j])

	for k in range(len(index)):
		file_name = index[k]
		fMRI_file = os.path.join(func_directory,file_name)
		os.system('mkdir '+output+'tmp/func')

        #Get slice and TR parameters from header
		info = os.popen('fslinfo '+fMRI_file)
		info = info.read()
		Slice = int(info[info.index('dim4')+6:info.index('datatype')-1])
		TR = float(info[info.index('pixdim4')+9:info.index('cal_max')-1])

		#Slice Timing
		if os.path.exists(output+'tmp/func/SliceTiming_'+file_name) == 0:
			if SliceOrder != 'NoSliceTiming':
				os.system('slicetimer -i '+fMRI_file+' -r '+str(TR)+' --ocustom='+SliceOrder+' -o '+output+'tmp/func/SliceTiming_'+file_name)
			else:
				os.system('cp '+fMRI_file+' '+output+'tmp/func/SliceTiming_'+file_name)
				
		#Remove the first N time points
		if os.path.exists(output+'tmp/func/rmftp_'+file_name) == 0:
			if first_TP == 0:
				tsize = Slice
				os.system('mv '+output+'tmp/func/SliceTiming_'+file_name+' '+output+'tmp/func/rmftp_'+file_name)
			if first_TP != 0:
				tsize = Slice - first_TP
				os.system('fslroi '+output+'tmp/func/SliceTiming_'+file_name+' '+output+'tmp/func/rmftp_'+file_name+' '+str(first_TP)+' '+str(tsize))

        #Motion correction
		if mc_method != 'None':
			if os.path.exists(output+'tmp/func/mc_'+file_name) == 0:
				Motion_corr(output+'tmp/func/rmftp_'+file_name,output,file_name,mc_method)
		else:
			os.system('mv '+output+'tmp/func/rmftp_'+file_name+' '+output+'tmp/func/mc_'+file_name)

        #Distortion field correction
		if fieldmap and fieldmap_method != '':
			tmp_name = Fieldmap_corr(output+'tmp/func/mc_'+file_name,output,data_directory,file_name,fieldmap_method)
		else:
			tmp_name = 'mc_'+file_name

        #Normalization
		if os.path.exists(output+'tmp/func/f2std_'+file_name) == 0:
			Normalization(output+'tmp/func/'+tmp_name,output,file_name,Standard_image_func,tsize,TR)

        #Smooth
		if fwhm != 0 and os.path.exists(output+'tmp/func/smooth_'+file_name)==0:
			tmp_name = Smooth(output+'tmp/func/f2std_'+file_name,output,file_name,fwhm)
		else:
			tmp_name = 'f2std_'+file_name

		os.system('mv '+output+'tmp/func/'+tmp_name+' '+output+'func/predenoise_'+file_name)

        #Denoise
		prefilter = []
		if denoise_method == 'no':
			prefilter.append('predenoise_'+file_name)

		if denoise_method == 'regress':
			if despike:
				os.system('matlab -nodesktop -nosplash -r "input=\''+output+'func/predenoise_'+file_name+'\'; toolbox_dir=\''+toolbox_dir+'\'; cd(toolbox_dir); x_Despike(input,toolbox_dir); exit"')
				os.system('matlab -nodesktop -nosplash -r "toolbox_dir=\''+toolbox_dir+'\'; cd(toolbox_dir); funcfile=\''+output+'tmp/func/despiked_wds.nii.gz\'; motionfile=\''+output+'func/mc_'+file_name+'.txt\'; motion=' \
					+str(motion)+'; CSFmask=\''+output+'Seg/'+sub+'_CSF_mask_'+Res+'mm.nii.gz\'; CSF='+str(CSF)+'; WMmask=\''+output+'Seg/'+sub+'_WM_mask_'+Res+'mm.nii.gz\'; WM='+str(WM)+ \
					'; x_extract_nuisance_covariates(funcfile,motionfile,motion,CSFmask,CSF,WMmask,WM); exit"')
				os.system('1dcat '+output+'tmp/func/despiked_wds_regress_out.txt > '+output+'func/despiked_'+file_name.replace('.nii.gz','_regress_out.1D'))
				os.system('3dTproject -ort '+output+'func/despiked_'+file_name.replace('.nii.gz','_regress_out.1D')+' -prefix '+output+'func/Regress_denoised_'+file_name+' -input '+output+'tmp/func/despiked_wds.nii.gz')
			else:
				os.system('matlab -nodesktop -nosplash -r "toolbox_dir=\''+toolbox_dir+'\'; cd(toolbox_dir); funcfile=\''+output+'func/predenoise_'+file_name+'\'; motionfile=\''+output+'func/mc_'+file_name+'.txt\'; motion=' \
					+str(motion)+'; CSFmask=\''+output+'Seg/'+sub+'_CSF_mask_'+Res+'mm.nii.gz\'; CSF='+str(CSF)+'; WMmask=\''+output+'Seg/'+sub+'_WM_mask_'+Res+'mm.nii.gz\'; WM='+str(WM)+ \
					'; x_extract_nuisance_covariates(funcfile,motionfile,motion,CSFmask,CSF,WMmask,WM); exit"')
				regress_txt = file_name.replace('.nii.gz','_regress_out.txt')
				regress_1D = regress_txt.replace('.txt','.1D')
				os.system('1dcat '+output+'func/predenoise_'+regress_txt+' > '+output+'func/'+regress_1D)
				os.system('3dTproject -ort '+output+'func/'+regress_1D+' -prefix '+output+'func/Regress_denoised_'+file_name+' -input '+output+'func/predenoise_'+file_name)

			prefilter.append('Regress_denoised_'+file_name)

		if denoise_method == 'ICA' and os.path.exists(output+'func/ICA_denoised_'+file_name) == 0:
			os.system('rm -r '+output+'tmp/func/ICA')
			os.system('python '+toolbox_dir+'ICA-AROMA/ICA_AROMA.py -i '+output+'func/predenoise_'+file_name+' -o '\
				+output+'tmp/func/ICA -mc '+output+'func/mc_'+file_name+'.txt -tr '+str(TR))
			os.system('mv '+output+'tmp/func/ICA/denoised_func_data_nonaggr.nii.gz '+output+'func/ICA_denoised_'+file_name)
			prefilter.append('ICA_denoised_'+file_name)

        #Filter or Detrend
		for l in range(len(prefilter)):
			if lowpass == 0 and highpass == 0:
				os.system('3dDetrend -polort 1 -prefix '+output+'tmp/func/filter_'+prefilter[l]+ \
					' '+output+'func/'+prefilter[l])
			else:
				os.system('3dFourier -lowpass '+str(lowpass)+' -highpass '+str(highpass)+' -prefix '+output+'tmp/func/filter_'+prefilter[l]+ \
					' '+output+'func/'+prefilter[l])
			os.system('fslmaths '+output+'func/predenoise_'+file_name+' -Tmean -add '+output+'tmp/func/filter_'+prefilter[l]+' '+output+'func/filter_'+prefilter[l])
			if os.path.exists(output+'func/filter_'+prefilter[l]) == 1:
				check = 1
			else:
				check = 0

		if check:
			os.system('rm -r '+output+'tmp/func')

	prep_files = os.listdir(output+'func')
	check = []
	for m in range(len(prep_files)):
		if prep_files[m].find('filter_')!=-1:
			check.append(prep_files[m]) 
	if len(check) == len(index) or len(check) == 2*len(index):
		os.system('rm -r '+output+'tmp')

pool = Pool(processes = 5)

for i in range(num,num+n):
	sub = sublist[i]
	pool.apply_async(prep_main,(sub,))

pool.close()
pool.join()

print('=======================Congratulations! Your job is done!=======================')







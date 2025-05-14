
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 15:30:40 2019

@author: reiters
"""

import cv2
import numpy as np
import moving_least_squares
import os.path
import h5py
import matplotlib.pyplot as plt
import pandas as pd
import glob
import time

def norm(data):
        with np.errstate(invalid='ignore'):
            return (data - np.mean(data)) / np.std(data)

def get_patch_at(img, pos, patch_size):
        x, y = pos
        t = np.array([[1., 0., patch_size - x], \
                      [0., 1., patch_size - y]], 'float32')
        patch = cv2.warpAffine(img, t, \
                (int(2 * patch_size + .5), int(2 * patch_size + .5)))
        return norm(patch)

def warp_image(img, p, q, alpha, grid_size):
    #img is trg
    #p is src
    #q is trg
    
    identity_maps = np.dstack(np.meshgrid( \
            *tuple(np.arange(0, s + grid_size, grid_size) \
                  for s in img.shape[1::-1]))).astype('float32')
    coords = identity_maps.reshape((-1, 2))
    mapped_coords = moving_least_squares.similarity( \
            p, q, coords, alpha=alpha)
    maps = mapped_coords.reshape(identity_maps.shape)
    
    mapped_coords_inv = moving_least_squares.similarity( \
            q, p, coords, alpha=alpha)
    maps_inv = mapped_coords_inv.reshape(identity_maps.shape)
    
    t = np.array([[grid_size, 0, 0], [0, grid_size, 0]], 'float32')
    maps = cv2.warpAffine(maps, t, img.shape[1::-1])
    maps_inv = cv2.warpAffine(maps_inv, t, img.shape[1::-1])
    
   
    return cv2.remap(img, maps, None, interpolation=cv2.INTER_LINEAR), maps


#%% Get the list of images to stitch and their initial matchfeatures files

#Put the 2 images in a folder, this is the folder address
dirname='/home/sam/pophale2023/example_data/'

#these are the 2 image names
file1='wakeImg.png' 
file2='sleepImg.png' 



src_mf=cv2.imread(dirname + '/' + file1,0)
trg_mf=cv2.imread(dirname + '/' + file2,0)

#%% Load up manual point correspondences and warp images
# to do full manual nonlinear matching, see stitching/manualStitch


mls_alpha=3.0
feat_file=h5py.File(dirname + '/' + file2 + '_matchFeatures','r')
tempPtsSrc=feat_file['tempPtsSrc'][:].astype('float32')
tempPtsTrg=feat_file['tempPtsTrg'][:].astype('float32')
feat_file.close()

trg_mf_warped,maps = warp_image(trg_mf, \
        tempPtsSrc, tempPtsTrg, mls_alpha, 32)

    

C = np.uint8(np.dstack((src_mf,trg_mf_warped,src_mf)))
plt.figure()
plt.imshow(C)


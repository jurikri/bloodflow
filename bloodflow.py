# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 13:17:24 2022

@author: PC
"""

#%%

from skimage import io # pip install scikit-image (체크필요)
import matplotlib.pyplot as plt
import numpy as np
# import cv2
from PIL import Image 
from tqdm import tqdm
from skimage.data import shepp_logan_phantom
from skimage.transform import radon, rescale
import os
import sys; 
sys.path.append('D:\\mscore\\code_lab\\')
sys.path.append('C:\\mscode')
import msFunction
import pickle

# filepath = 'C:\\Temp\\Image5_ time series_ spon_line scan_Galvano.tif'

filepath = 'D:\\HR_bloodflow\\STZ\\3C1\\' + '3C1_Image2_line scan.bmp'

filepath = 'D:\\HR_bloodflow\\STZ\\3C1\\' + '3C1_Image3_line scan.bmp'

filepath = 'D:\\HR_bloodflow\\STZ\\3C1\\' + '3C1_Image7_line scan.bmp'


#%% auto filepath

path = 'D:\\HR_bloodflow\\'
total_path = []
for (path, dir, files) in os.walk("D:/"):
    for filename in files:
        ext = os.path.splitext(filename)[-1]
        if ext == '.bmp':
            print("%s/%s" % (path, filename), path[16:19])
            
            filepath = path + '//' + filename
            total_path.append([filepath, path[16:19]])

#%%        

mssave = msFunction.msarray([len(total_path)])

for i in tqdm(range(0, len(total_path))):

    filepath = total_path[i][0]
    img0 = io.imread(filepath) # 3 dimensions : frames x width x height
    if img0.shape[0] < img0.shape[1]: img0 = np.transpose(img0)
    
    # plt.imshow(img0)
    
    if False:
        for j in range(5):
            plt.figure(); plt.imshow(img0[int(img0.shape[1]*(7*j)):int(img0.shape[1]*(7*(j+1))), :])
            
    width = img0.shape[1]
    bins = int(width/2)
    msbins = np.arange(0, img0.shape[0]-bins+1, bins, dtype=int)
    peak_save = []
    
    method1 = False
    method2 = True
    
    for row in msbins:
        stdsave, stdsave2 = [], []
        resolution = 1
        angle_list = np.arange(0, 180.0001, resolution)
        figsw = False
        for ro in angle_list:
            crop = img0[row :row + width, :]
            crop_img = Image.fromarray(crop / np.mean(crop, axis=0) - 1)
            rotate_img = crop_img.rotate(ro)
            rotate_img_array = np.array(rotate_img, dtype=float)
            
            crop2 = np.array(rotate_img_array)
            
            if np.isnan(np.mean(crop2)): print('nan')
            
            # method1
            if method1:
                theta = np.linspace(0., 180., max(crop2.shape), endpoint=False)
                sinogram = radon(crop2, theta=theta)
                score = np.std(sinogram) / np.std(crop2)
                stdsave.append(score)
            
            # method2
            
            mean_trace = np.nanmean(crop2, axis=0)
            formula2 = np.std(mean_trace)
            stdsave2.append(formula2)
            
        if method1:
            stdsave = np.array(stdsave)
            ws = int(round(31 * (0.1 / resolution)))
            smooth = np.convolve(stdsave, np.ones((ws,))/ws, mode='valid')
            stdsave_smooth = np.zeros(stdsave.shape)
            stdsave_smooth[int(ws/2):int(ws/2)+len(smooth)] = smooth
            plt.figure()
            plt.plot(stdsave[1:-1])
            plt.plot(stdsave_smooth[1:-1])
            
        if method2:
            stdsave2 = np.array(stdsave2)
            ws = int(round(31 * (0.1 / resolution)))
            smooth = np.convolve(stdsave2, np.ones((ws,))/ws, mode='valid')
            stdsave_smooth = np.zeros(stdsave2.shape)
            stdsave_smooth[int(ws/2):int(ws/2)+len(smooth)] = smooth
            
        mix = np.argmax(stdsave_smooth)
        ms_angle = angle_list[mix]
        ms_value = stdsave2[mix]
        
        peak_save.append([ms_angle, ms_value])
    peak_save = np.array(peak_save)

    median_angle = np.median(np.abs(peak_save[:,0] - 90))
    print(median_angle)
    msid = os.path.basename(filepath)[:-4]
    mssave[i] = [peak_save, msid]
    
    if False:
        plt.figure();
        plt.title(str(i) + '_' +  filepath)
        plt.scatter(list(range(len(peak_save[:,0]))), np.abs(peak_save[:,0] - 90), s=5)
        plt.plot(np.ones(len(peak_save[:,0])) * median_angle, c='r')

#%%
psave = 'D:\\HR_bloodflow\\blood_flow_save.pickle'
if not(os.path.isfile(psave)) or True:
    with open(psave, 'wb') as f:  # Python 3: open(..., 'rb')
        pickle.dump(mssave, f, pickle.HIGHEST_PROTOCOL)
        print(psave, '저장되었습니다.')
    

#%% vis
psave = 'D:\\HR_bloodflow\\blood_flow_save.pickle'
with open(psave, 'rb') as file:
    mssave = pickle.load(file)

for i in range(len(mssave)):
    filepath = total_path[i][0]
    peak_save = mssave[i][0]
    
    median_angle = np.median(np.abs(peak_save[:,0] - 90))
    plt.figure();
    plt.title(str(i) + '_' +  filepath)
    plt.scatter(list(range(len(peak_save[:,0]))), np.abs(peak_save[:,0] - 90), s=5)
    plt.plot(np.ones(len(peak_save[:,0])) * median_angle, c='r')

#%% grouping
import pandas as pd
filepath = 'C:\\SynologyDrive\\worik in progress\\20220517 - flow calc\\' + 'STZ 정리_ms.xlsx'

groupinfo = []

df = np.array(pd.read_excel(filepath, sheet_name=0, names=None))
for i in range(len(df)):
    msid = df[i,0][:-4]
    quality = df[i,2]
    
    size = df[i,7]
    spatial = float(size[:size.find(' [um]')]) # um
    time = float(size[size.find('*')+2 : size.find(' [ms]')]) # um
    label = 'STZ'
    se = df[i,9]
    SE = df[i,10]
    
    groupinfo.append([msid, quality, spatial, time, label, SE, se])
    
df = np.array(pd.read_excel(filepath, sheet_name=1, names=None))
for i in range(len(df)):
    msid = df[i,0][:-4]
    quality = df[i,2]
    
    size = df[i,7]
    spatial = float(size[:size.find(' [um]')]) # um
    time = float(size[size.find('*')+2 : size.find(' [ms]')]) # um
    label = 'Vehicle'
    se = df[i,9]
    SE = df[i,10]
    
    groupinfo.append([msid, quality, spatial, time, label, SE, se])


#%% post analysis
idindex = np.array(mssave)[:,1]

mssave7 = []
for i in range(len(groupinfo)):

    msid = groupinfo[i][0] 
    idix = np.where(idindex == msid)[0][0] 
    msplot = mssave[idix][0][:,0]
    mean_speed = np.mean(msplot)
    mssave7.append(groupinfo[i] + [mean_speed])
        
    # plt.scatter(list(range(len(msplot))), msplot, s=2)
mssave7 = np.array(mssave7)

mssave8 = msFunction.msarray([2])
SElist = list(set(np.array(mssave7[:,5], dtype=float)))
for SE in SElist:
    SEix = np.logical_and(np.array(mssave7[:,1], dtype=float) <= 1, np.array(mssave7[:,5], dtype=float) == SE)
    if len(np.where(SEix)[0]) == 2:
        before_ix = np.where(np.logical_and(SEix, np.array(mssave7[:,6], dtype=float) == 1))[0][0]
        after_ix = np.where(np.logical_and(SEix, np.array(mssave7[:,6], dtype=float) == 2))[0][0]

        before_speed = float(mssave7[before_ix,7]) * (float(mssave7[before_ix,2]) / float(mssave7[before_ix,3]))
        after_speed = float(mssave7[after_ix,7]) * (float(mssave7[after_ix,2]) / float(mssave7[after_ix,3]))

        print(before_ix, after_ix)
        
        label = None
        if mssave7[before_ix,4] == 'Vehicle': label = 0
        if mssave7[before_ix,4] == 'STZ': label = 1
        
        mssave8[label].append([before_speed, after_speed])
        
Aprism1 = np.array(mssave8[0])
Aprism2 = np.array(mssave8[1])

#%%


# for i in range(len(mix)):
    
for i in range(0, 50):
    row = msbins[i]
    otimized_angle = np.round(angle_list[int(peak_save[i,0])], 2)
    
    crop = img0[row:row+width, :]
    
    # crop = img0[1230:1500, :]
    # plt.imshow(crop)
    
    line_matrix = np.zeros(crop.shape)
    lix = np.arange(10, crop.shape[1], 10, dtype=int)
    line_matrix[:,lix] = np.max(crop) * 2
    line_matrix = Image.fromarray(line_matrix)
    
    rotate_img = line_matrix.rotate(-otimized_angle)
    rotate_rgb = rotate_img.convert("RGB")
    
    crop_img = Image.fromarray(crop * 4)
    crop_rgb = crop_img.convert("RGB")
    crop_array = np.array(crop_rgb)
    
    crop_array[:,:,0] = np.array(rotate_rgb)[:,:,0]
    crop_array_img =  Image.fromarray(crop_array)
    
    fig = plt.figure()
    rows = 1
    cols = 2
    
    ax1 = fig.add_subplot(rows, cols, 2)
    ax1.imshow(crop_array_img)
    ax1.set_title('Merge ' + str(otimized_angle) + ' degree')
    ax1.axis("off")
    
    ax2 = fig.add_subplot(rows, cols, 1)
    ax2.imshow(Image.fromarray(crop))
    ax2.set_title('Original_# ' + str(i))
    ax2.axis("off")
    


#%% radon transform


import numpy as np
import matplotlib.pyplot as plt
from skimage.data import shepp_logan_phantom
from skimage.transform import radon, rescale


theta = np.linspace(0., 180., max(crop.shape), endpoint=False)
sinogram = radon(crop, theta=theta)

plt.figure()
plt.imshow(crop)

plt.figure()
plt.imshow(sinogram)


resolution = 1
angle_list = np.arange(0, 180.0001, resolution)
figsw = False

stdsave, stdsave2 = [], []
for ro in angle_list:
    crop_img = Image.fromarray(crop - np.mean(crop, axis=0))
    rotate_img = crop_img.rotate(ro)
    rotate_img_array = np.array(rotate_img, dtype=float)
    rotate_img_array[rotate_img_array==0] = np.nan

    if figsw:
        plt.figure()
        plt.imshow(np.array(rotate_img))
        plt.title(str(ro))
        
    vix = np.where(np.isnan(np.mean(rotate_img_array, axis=0))==0)[0]
    theta = np.linspace(0., 180., max(rotate_img_array.shape), endpoint=False)
    sinogram = radon(rotate_img_array, theta=theta)

    score = np.nanstd(rotate_img_array) / np.nanstd(sinogram)
    stdsave.append(score)
    
    mean_trace = np.nanmean(rotate_img_array[:,vix], axis=0)
    formula2 = np.std(mean_trace)
    stdsave2.append(formula2)
    
    #%%
    
otimized_angle = 159

otimized_angle = 169

line_matrix = np.zeros(crop.shape)
lix = np.arange(10, crop.shape[1], 10, dtype=int)
line_matrix[:,lix] = np.max(crop) * 2
line_matrix = Image.fromarray(line_matrix)

rotate_img = line_matrix.rotate(-otimized_angle)
rotate_rgb = rotate_img.convert("RGB")

crop_img = Image.fromarray(crop * 4)
crop_rgb = crop_img.convert("RGB")
crop_array = np.array(crop_rgb)

crop_array[:,:,0] = np.array(rotate_rgb)[:,:,0]
crop_array_img =  Image.fromarray(crop_array)

fig = plt.figure()
rows = 1
cols = 2

ax1 = fig.add_subplot(rows, cols, 2)
ax1.imshow(crop_array_img)
ax1.set_title('Merge ' + str(otimized_angle) + ' degree')
ax1.axis("off")

ax2 = fig.add_subplot(rows, cols, 1)
ax2.imshow(Image.fromarray(crop))
ax2.set_title('Original_# ' + str(i))
ax2.axis("off")


#%% toy sample - test1


def ms_toy_img_gen(angle=None, width=30, line_resolution=3):
    import numpy as np
    from PIL import Image 
    hw = int(width / np.sqrt(2) / 2)
    line_matrix = np.zeros((width*3, width))
    lix = np.arange(0, line_matrix.shape[1], line_resolution, dtype=int)
    line_matrix[:,lix] = 10
    line_matrix = Image.fromarray(line_matrix)
    rotate_img = line_matrix.rotate(-angle)
    rotate_rgb = rotate_img.convert("RGB")
    crop = np.array(rotate_rgb)[:,:,0]
    crop = crop[int(width*3/2) - hw : int(width*3/2) + hw, int(width/2) - hw : int(width/2) + hw]
    
    return crop

def ms_randon_transform(img=None):
    ar = 10
    angle_list = np.arange(30, 200.0001, ar)
    stdsave = []
    crop = np.array(img)
    for ro in angle_list:
        
        crop_img = Image.fromarray(crop / np.mean(crop, axis=0) - 1)
        rotate_img = crop_img.rotate(ro)
        rotate_img_array = np.array(rotate_img, dtype=float)
        crop2 = np.array(rotate_img_array)
        
        theta = np.linspace(0., 180., max(crop2.shape), endpoint=False)
        sinogram = radon(crop2, theta=theta)
        score = np.nanstd(sinogram) / np.nanstd(crop2)
        stdsave.append(score)
    
    return angle_list, stdsave

def ms_ng_angle_method(img=None):
    import numpy as np
    
    ar = 10
    angle_list = np.arange(30, 200.0001, ar)
    stdsave = []
    crop = np.array(img)
    for ro in angle_list:
        
        crop_img = Image.fromarray(crop / np.mean(crop, axis=0) - 1)
        rotate_img = crop_img.rotate(ro)
        rotate_img_array = np.array(rotate_img, dtype=float)
        crop2 = np.array(rotate_img_array)
        
        mean_trace = np.nanmean(crop2, axis=0)
        formula2 = np.std(mean_trace)
        stdsave.append(formula2)
        
    stdsave = np.array(stdsave)
    ws = int(round(31 * (0.1 / ar)))
    if ws > 0:
        smooth = np.convolve(stdsave2, np.ones((ws,))/ws, mode='valid')
        stdsave_smooth = np.zeros(stdsave2.shape)
        stdsave_smooth[int(ws/2):int(ws/2)+len(smooth)] = smooth
    else: stdsave_smooth = np.array(stdsave)
    
    return angle_list, stdsave_smooth

def ms_minmax(X): # [0,1]
    X = np.array(X)
    msmin = np.min(X)
    msmax = np.max(X)
    
    X_std = (X - msmin) / (msmax - msmin)
    X_scaled = X_std * (1 - 0) + 0
    return X_scaled


time_resolution = 10
angle_list = np.arange(30, 200.0001, time_resolution)

ro = 40
img = ms_toy_img_gen(angle=ro)
plt.imshow(img)

xaxis, stdsave = ms_randon_transform(img=img)
angle_list, stdsave_smooth = ms_ng_angle_method(img=img)

plt.figure()
plt.plot(xaxis, ms_minmax(stdsave_smooth), label = 'ngangle')
plt.plot(xaxis, ms_minmax(stdsave), label = 'randon')
plt.legend()

# test2
# test_img = img0[int(img0.shape[1]*(7*j)):int(img0.shape[1]*(7*(j+1))), :]
test_img = np.load('C:\\SynologyDrive\\worik in progress\\20220517 - flow calc\\test_img.npy')        
plt.imshow(test_img)

n = 5
w = test_img.shape[1]
test_img2 = test_img[w*n:w*(n+1),:]
plt.imshow(test_img2)

img = np.array(test_img2)
xaxis, stdsave = ms_randon_transform(img=img)
angle_list, stdsave_smooth = ms_ng_angle_method(img=img)

plt.figure()
plt.plot(xaxis, ms_minmax(stdsave_smooth), label = 'ngangle')
plt.plot(xaxis, ms_minmax(stdsave), label = 'randon')
plt.legend()

# test3
filepath = 'C:\\SynologyDrive\\worik in progress\\20220517 - flow calc\\2022-03-24_TR_VNS_5m-line.tif'
img0 = io.imread(filepath) # 3 dimensions : frames x width x height
test_img = np.array(img0)[:,:,0][:,120:330]

n = 5
w = test_img.shape[1]
test_img2 = test_img[w*n:w*(n+1),:]
plt.imshow(test_img2)

img = np.array(test_img2)
xaxis, stdsave = ms_randon_transform(img=img)
angle_list, stdsave_smooth = ms_ng_angle_method(img=img)

plt.figure()
plt.plot(xaxis, ms_minmax(stdsave_smooth), label = 'ngangle')
plt.plot(xaxis, ms_minmax(stdsave), label = 'randon')
plt.legend()














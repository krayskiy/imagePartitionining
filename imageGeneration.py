import matplotlib.pyplot as plt
import numpy as np
import random
import skimage

# returns alteredImg, mask

def add_noise(image,prob):
    '''
    Add salt and pepper noise to image
    prob: Probability of the noise
    '''
    output = np.zeros(image.shape,np.uint8)
    thres = 1 - prob 
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            rdn = random.random()
            if rdn < prob:
                output[i][j] = 0
            elif rdn > thres:
                output[i][j] = 255
            else:
                output[i][j] = image[i][j]
    return output

def genAndGetMask(img, width, height):
    hmask = np.full((29, 30), 255)
    vmask = np.full((30, 29), 255)

    numEdges = random.randint(2, 5)*2 #numShapes * 2
    xBorder = set()
    yBorder = set()
    while (len(xBorder) < numEdges):
        xBorder.add(random.randint(1, width-1))
    while (len(yBorder) < numEdges):
        yBorder.add(random.randint(1, height-1))
    xBorder = sorted(xBorder)
    yBorder = sorted(yBorder)
    print("xBorder", xBorder)
    print("yBorder", yBorder)
    idx = 0
    while (idx < len(xBorder)):
        # add shapes to img
        img[xBorder[idx]:xBorder[idx+1]+1, yBorder[idx]:yBorder[idx+1]+1] = 0

        # draw edges of mask
        # mask[xBorder[idx]:xBorder[idx+1]+1, yBorder[idx]] = 0
        # mask[xBorder[idx]:xBorder[idx+1]+1, yBorder[idx+1]] = 0
        # mask[xBorder[idx], yBorder[idx]:yBorder[idx+1]+1] = 0
        # mask[xBorder[idx+1], yBorder[idx]:yBorder[idx+1]+1] = 0

        #left edge of mask
        vmask[xBorder[idx]-1, yBorder[idx]-1:yBorder[idx+1]] = 0
        #right edge of mask
        vmask[xBorder[idx+1], yBorder[idx]-1:yBorder[idx+1]] = 0
        #top edge of mask
        hmask[xBorder[idx]-1:xBorder[idx+1], yBorder[idx]-1] = 0
        #bottom edge of mask
        hmask[xBorder[idx]-1:xBorder[idx+1], yBorder[idx+1]] = 0

        idx += 2

    return img, vmask, hmask


def generateAndSave(name: str):
    piecewise = np.full((30, 30), 255)
    piecewise, vmask,  hmask = genAndGetMask(piecewise, 30, 30)
    #generate noise
    piecewise = add_noise(piecewise, 0.05)
    pname = "images/"+ name + ".png"
    vmaskname = "images/"+ name + "vmask.png"
    hmaskname = "images/"+ name + "hmask.png"
    plt.imsave(pname, piecewise, cmap='gray', vmin=0, vmax=255)
    plt.imsave(vmaskname, vmask, cmap='gray', vmin=0, vmax=255)
    plt.imsave(hmaskname, hmask, cmap='gray', vmin=0, vmax=255)
    # plt.imshow(piecewise, cmap='gray', vmin=0, vmax=255)
    plt.imshow(piecewise, cmap='gray', vmin=0, vmax=255)
    plt.imshow(vmask, cmap='gray', vmin=0, vmax=255)
    plt.imshow(hmask, cmap='gray', vmin=0, vmax=255)
    # piecewise = np.ones((30, 30), dtype=np.uint8)



for i in range(1,11):
    generateAndSave(str(i))
import matplotlib.pyplot as plt
import numpy as np
import random
import skimage

# returns alteredImg, mask


def alterAndGetMask(img, width, height):
    mask = img.copy()
    numEdges = random.randint(2, 4)*2
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
        mask[xBorder[idx]:xBorder[idx+1]+1, yBorder[idx]] = 0
        mask[xBorder[idx]:xBorder[idx+1]+1, yBorder[idx+1]] = 0
        mask[xBorder[idx], yBorder[idx]:yBorder[idx+1]+1] = 0
        mask[xBorder[idx+1], yBorder[idx]:yBorder[idx+1]+1] = 0

        idx += 2
    return img, mask


# def plotnoise(img, mode, r, c, i):
#     plt.subplot(r, c, i)
#     if mode is not None:
#         gimg =
#         plt.imshow(gimg)
#     else:
#         plt.imshow(img)
#     plt.title(mode)
#     plt.axis("off")


piecewise = np.full((30, 30), 255)
# piecewise = np.zeros((30, 30), dtype=np.uint8)
piecewise, mask = alterAndGetMask(piecewise, 30, 30)
# piecewise = skimage.util.random_noise(piecewise, mode="gaussian")
noise = (np.random.standard_normal((30, 30))*256 /
         6*(1 if random.random() < 0.5 else -1))+128
print(noise)
piecewise += noise.astype(np.uint8)

# piecewise += noise.astype(np.uint8)
plt.imsave('1.png', piecewise, cmap='gray', vmin=0, vmax=255)
plt.imsave('2.png', mask, cmap='gray', vmin=0, vmax=255)
# plt.imshow(piecewise, cmap='gray', vmin=0, vmax=255)

# piecewise = np.ones((30, 30), dtype=np.uint8)


plt.imshow(piecewise, cmap='gray', vmin=0, vmax=255)

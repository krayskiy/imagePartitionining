import matplotlib.pyplot as plt
import numpy as np
import time


def generate_piecewise(size: int, shapes: int = 3):
    rng = np.random.default_rng()
    start = 0
    colors = []
    regions = shapes + 1
    for r in range(regions):
        start += rng.integers(128//regions, 255//regions)
        colors.append(start)
    rng.shuffle(colors)
    piecewise = np.ones((size, size))
    piecewise *= colors[0]
    for i in range(1, regions):
        start_x = rng.uniform(0, size-3)
        start_y = rng.uniform(0, size-3)
        end_x = start_x + rng.normal((size-start_x)/2)
        end_x = max(end_x, start_x+3)
        end_y = start_y + rng.normal((size-start_y)/2)
        end_y = max(end_y, start_y+3)
        start_x, start_y, end_x, end_y = np.array([start_x, start_y, end_x, end_y],dtype = int)
        piecewise[start_y:end_y, start_x:end_x] = colors[i]
    masks = []
    for i in range(1, regions):
        masks.append(piecewise==colors[i])
    return piecewise, masks


def add_noise(piecewise: np.array, intensity: int = 8):
    rng = np.random.default_rng()
    noise = rng.normal(scale=intensity, size=piecewise.shape)
    piecewise += noise.astype(int)
    return np.clip(piecewise, 0, 255)

def get_edge_masks(masks):
    imheight, imwidth = masks[0].shape
    hmask = np.zeros((imheight - 1, imwidth))
    vmask = np.zeros((imheight, imwidth-1))
    for mask in masks:
        mask = mask.astype(int)
        for i in range(imheight):
            for j in range(imwidth):
                try:
                    hdif = mask[i+1, j] - mask[i,j]
                    vdif = mask[i,j+1] - mask[i,j]
                    if hdif != 0:
                        hmask[i,j] = 1
                    if vdif != 0:
                        vmask[i,j] = 1
                except IndexError:
                    continue 
    return hmask, vmask


piecewise, masks = generate_piecewise(30)
hmask, vmask = get_edge_masks(masks)
piecewise = add_noise(piecewise)
current = int(time.time())

# print("hmask shape = ", hmask.shape)
# print("vmask shape = ", vmask.shape)
# fig,(ax1,ax2,ax3) = plt.subplots(1,3)
# ax1.imshow(piecewise, cmap='gray', vmin=0, vmax=255)
# ax2.imshow(hmask, cmap='gray', vmin=0, vmax=1)
# ax3.imshow(vmask, cmap='gray', vmin=0, vmax=1)
# fig.show()


for i in range(1,11):
    piecewise, masks = generate_piecewise(30)
    hmask, vmask = get_edge_masks(masks)
    piecewise = add_noise(piecewise)
    plt.imsave(f"./images2/im{i}.png", piecewise, cmap='gray', vmin=0, vmax=255)
    plt.imsave(f"./images2/im{i}_hmask.png", hmask, cmap='gray', vmin=0, vmax=1)
    plt.imsave(f"./images2/im{i}_vmask.png", vmask, cmap='gray', vmin=0, vmax=1)



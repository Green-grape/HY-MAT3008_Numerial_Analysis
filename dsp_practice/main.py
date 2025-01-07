import os

import cv2
import numpy as np


def convolution2d(image, kernel, kernel_size=3):
    if kernel.shape != (kernel_size, kernel_size):
        raise ValueError("Kernel size must be square")

    h, w = image.shape

    pad_size = kernel_size // 2
    padded_image = np.pad(image, pad_size)
    output = np.zeros((h, w), dtype=np.float64)

    for i in range(h):
        for j in range(w):
            output[i, j] = np.sum(
                padded_image[i : i + kernel_size, j : j + kernel_size] * kernel
            )
    return output


def gaussian_smooth(image, kernel_size=3, sigma=1.0):
    kernel = cv2.getGaussianKernel(kernel_size, sigma)
    kernel = kernel @ kernel.T
    return convolution2d(image, kernel, kernel_size)


def sobel_filter(image):
    sobel_x = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]], dtype=np.float64)
    sobel_y = np.array([[-1, -2, -1], [0, 0, 0], [1, 2, 1]], dtype=np.float64)

    sobel_x_filtered = convolution2d(image, sobel_x)
    sobel_y_filtered = convolution2d(image, sobel_y)
    sobel_combined = np.sqrt(sobel_x_filtered**2 + sobel_y_filtered**2)

    return sobel_combined.astype(np.uint8)


if __name__ == "__main__":
    cur_dir_path = os.path.dirname(os.path.realpath(__file__))
    img = cv2.imread(f"{cur_dir_path}/image2.jpg", cv2.IMREAD_GRAYSCALE)
    gs_img = gaussian_smooth(img).astype(np.uint8)
    sobel_img = sobel_filter(img).astype(np.uint8)

    custom_kernel = np.array([[0, -1, 0], [-1, 5, -1], [0, -1, 0]])  # Sharpening kernel
    sharpened_img = convolution2d(img, custom_kernel).astype(np.uint8)

    result = cv2.hconcat([img, gs_img, sobel_img, sharpened_img])

    cv2.imshow("Result", result)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

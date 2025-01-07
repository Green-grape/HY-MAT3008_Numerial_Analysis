import os
from typing import List, Tuple

import cv2
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def compute_residuals(
    params: List[int], pts1: List[Tuple[float, float]], pts2: List[Tuple[float, float]]
):
    residuals = []
    for p1, p2 in zip(pts1, pts2):
        x1, y1 = p1
        x2, y2 = p2
        norm = params[6] * x1 + params[7] * y1 + 1
        x1_proj = (params[0] * x1 + params[1] * y1 + params[2]) / norm
        y1_proj = (params[3] * x1 + params[4] * y1 + params[5]) / norm
        residuals.append((x2 - x1_proj) ** 2)
        residuals.append((y2 - y1_proj) ** 2)
    return np.array(residuals)


def compute_jacobian(params, pts1, pts2):
    jacobian = []
    for p1, p2 in zip(pts1, pts2):
        x1, y1 = p1
        x2, y2 = p2
        norm = params[6] * x1 + params[7] * y1 + 1
        x1_proj = (params[0] * x1 + params[1] * y1 + params[2]) / norm
        y1_proj = (params[3] * x1 + params[4] * y1 + params[5]) / norm
        d_x_proj = [
            2 * (x1_proj - x2) * x1 / norm,
            2 * (x1_proj - x2) * y1 / norm,
            2 * (x1_proj - x2) / norm,
            0,
            0,
            0,
            -2 * (x1_proj - x2) * x1_proj * x1 / norm,
            -2 * (x1_proj - x2) * x1_proj * y1 / norm,
        ]
        d_y_proj = [
            0,
            0,
            0,
            2 * (y1_proj - y2) * x1 / norm,
            2 * (y1_proj - y2) * y1 / norm,
            2 * (y1_proj - y2) / norm,
            -2 * (y1_proj - y2) * y1_proj * x1 / norm,
            -2 * (y1_proj - y2) * y1_proj * y1 / norm,
        ]
        jacobian.append(d_x_proj)
        jacobian.append(d_y_proj)
    return np.array(jacobian)


def compute_homography(pts1, pts2):
    params = np.array([1, 0, 0, 0, 1, 0, 0, 0], dtype=np.float64)

    max_iter = 1000
    lambda_ = 1e-3
    tolerance = 1e-3

    for i in range(max_iter):
        residuals = compute_residuals(params, pts1, pts2)
        jacobian = compute_jacobian(params, pts1, pts2)
        hessian = jacobian.T @ jacobian
        hessian_lm = hessian + lambda_ * np.eye(hessian.shape[0])

        gradient = jacobian.T @ residuals
        delta = np.linalg.solve(hessian_lm, -gradient)

        if np.linalg.norm(delta) < tolerance:
            print(f"Converged after {i} iterations")
            break

        residuals_new = compute_residuals(params + delta, pts1, pts2)
        if np.linalg.norm(residuals_new) < np.linalg.norm(residuals):
            params += delta
            lambda_ /= 5
        else:
            lambda_ *= 5
    H_optimized = np.concatenate((params, np.array([1])), axis=0).reshape(3, 3)
    return H_optimized


if __name__ == "__main__":
    cur_dir_path = os.path.dirname(os.path.realpath(__file__))
    img1 = cv2.imread(
        os.path.join(cur_dir_path, "./Figure_8.png"), cv2.IMREAD_GRAYSCALE
    )

    height, width = img1.shape
    src_points = np.float32(
        [[0, 0], [width - 1, 0], [0, height - 1], [width - 1, height - 1]]
    )
    dst_points = np.float32(
        [
            [0, height / 3],  # 왼쪽 위
            [width - 1, 0],  # 오른쪽 위
            [0, height * 2 / 3],  # 왼쪽 아래
            [width - 1, height - 1],  # 오른쪽 아래
        ]
    )
    H, _ = cv2.findHomography(src_points, dst_points)
    wraped_img = cv2.warpPerspective(img1, H, (width, height))

    cv2.imshow("Warped Image", wraped_img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

    cv2.imshow("Original Image", img1)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

    orb = cv2.ORB_create()
    kp1, des1 = orb.detectAndCompute(img1, None)
    kp2, des2 = orb.detectAndCompute(wraped_img, None)

    bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)
    matches = bf.match(des1, des2)
    matches = sorted(matches, key=lambda x: x.distance)

    print(f"Number of matches: {len(matches)}")
    cv2.imshow(
        "Matches",
        cv2.drawMatches(
            img1,
            kp1,
            wraped_img,
            kp2,
            matches[:10],
            None,
            flags=cv2.DrawMatchesFlags_NOT_DRAW_SINGLE_POINTS,
        ),
    )
    cv2.waitKey(0)
    cv2.destroyAllWindows()

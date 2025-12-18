import numpy as np
import matplotlib.pyplot as plt

def get_equal_spaced_ellipse_points(a, b, theta1, theta2, ds):
    def f(theta):
        # dtheta/ds
        return 1.0 / np.sqrt((a**2 * np.sin(theta)**2) + (b**2 * np.cos(theta)**2))

    points = []
    current_theta = theta1
    direction = 1 if theta2 > theta1 else -1
    step = ds * direction

    # Initial Point
    points.append((a * np.cos(current_theta), b * np.sin(current_theta)))

    while (direction * current_theta) < (direction * theta2):
        # RK4 Integration Step
        k1 = step * f(current_theta)
        k2 = step * f(current_theta + k1/2)
        k3 = step * f(current_theta + k2/2)
        k4 = step * f(current_theta + k3)

        current_theta += (k1 + 2*k2 + 2*k3 + k4) / 6

        if (direction * current_theta) > (direction * theta2):
            current_theta = theta2

        points.append((a * np.cos(current_theta), b * np.sin(current_theta)))
        if current_theta == theta2: break

    return np.array(points)

def plot_ellipse_points(points, a, b):
    plt.figure(figsize=(12, 6))

    # Plot the segments
    plt.plot(points[:, 0], points[:, 1], 'o-', markersize=4, label='Equidistant Points')

    # Highlight the first and last points
    plt.plot(points[0, 0], points[0, 1], 'ro', label='Start')
    plt.plot(points[-1, 0], points[-1, 1], 'go', label='End')

    plt.title(f"Elliptical Arc (a={a}, b={b}) with Equal Arc-Length Segments")
    plt.xlabel("X")
    plt.ylabel("Y")

    # CRITICAL: Set equal aspect ratio
    plt.gca().set_aspect('equal', adjustable='box')

    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.show()

# Parameters
A, B = 10.0, 1.0
T1, T2 = -np.pi/2, np.pi/2  # Half the ellipse
DS = 0.2           # Length of each segment

pts = get_equal_spaced_ellipse_points(A, B, T1, T2, DS)
plot_ellipse_points(pts, A, B)


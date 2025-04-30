import numpy as np

Lx = 15.0
Ly = 12.0
tolerance = 1e-4
for multiple in range(1, 10): # 网格分辨率增大且保持 dx = dy
    nx = 5 * multiple + 1
    ny = 4 * multiple + 1
    max_iter = 10000
    least_iter = 10000
    dx = Lx / (nx - 1)
    dy = Ly / (ny - 1)
    beta = dx / dy
    best_omega = 0

    # 修改omega
    start = 0
    end = 2
    step = 0.05
    for m in range(1, int((end - start) / step) + 1):
        omega = start + m * step
        if omega >= 2:
            break
        # 初始化（保持边界条件不变）
        T = np.full((nx, ny), 20.0)
        T[:, ny - 1] = 100.0

        # SOR迭代
        for iter in range(max_iter):
            max_diff = 0.0
            # 按行优先顺序遍历内部节点
            for i in range(1, nx - 1):
                for j in range(1, ny - 1):
                    old_val = T[i, j]
                    T[i, j] = T[i, j] + omega * 0.5 / (1 + beta * beta) * (
                                T[i + 1, j] + T[i - 1, j] + beta * beta * (T[i, j + 1] + T[i, j - 1]) - 2 * (
                                    1 + beta * beta) * T[i, j])
                    max_diff = max(max_diff, abs(T[i, j] - old_val))
            # 收敛判断
            if max_diff < tolerance:
                if iter < least_iter:
                    least_iter = iter
                    best_omega = omega
                    break

    # 输出每个(dx, dy)下最优omega
    print(f"dx={dx}，dy={dy}时，最优ω={best_omega}")

import numpy as np
import matplotlib.pyplot as plt

omega = 1.95  # 松弛因子

Lx = 15.0
Ly = 12.0
nx = 51
ny = 41
max_iter = 10000 # 最大迭代次数
tolerance = 1e-4
dx = Lx / (nx - 1)
dy = Ly / (ny - 1)
beta = dx / dy

# 初始化（保持边界条件不变）
T = np.full((nx, ny), 20.0)
T[ :, ny - 1] = 100.0

# SOR迭代
for iter in range(max_iter):
    max_diff = 0.0
    # 按行优先顺序遍历内部节点
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            old_val = T[i, j]
            T[i, j] = T[i, j] + omega * 0.5 / (1 + beta * beta) * (T[i+1, j] + T[i-1, j] + beta * beta * (T[i, j+1] + T[i, j-1]) - 2 * (1 + beta * beta) * T[i, j])
            max_diff = max(max_diff, abs(T[i, j] - old_val))
    # 收敛判断
    if max_diff < tolerance:
        print(f"SOR收敛于第{iter}次迭代，ω={omega}，最大变化{max_diff:.6f}")
        break
else:
    print("未收敛")

# 可视化代码
x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x, y)

plt.figure(figsize=(10, 6))
plt.contourf(X, Y, T.T, levels=20, cmap='jet')
plt.colorbar(label='T (°C)')
plt.title(f"SOR method temperature distribution (ω={omega})")
plt.xlabel("x (cm)")
plt.ylabel("y (cm)")
plt.show()
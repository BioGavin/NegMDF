import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt

# 示例点集
np.random.seed(42)
points = np.random.rand(10, 2)  # 生成10个二维点（随机）

# 计算凸包
hull = ConvexHull(points)

# 打印凸包顶点的索引
print("Convex hull vertices:", hull.vertices)

# 绘制点和凸包
plt.scatter(points[:, 0], points[:, 1], label="Points")  # 绘制所有点
for simplex in hull.simplices:  # 遍历凸包边
    plt.plot(points[simplex, 0], points[simplex, 1], 'r-')  # 绘制凸包边
plt.legend()
# 保存为 PNG 图片
plt.savefig("convex_hull_plot.png", dpi=300, bbox_inches="tight")  # 高分辨率保存
plt.show()
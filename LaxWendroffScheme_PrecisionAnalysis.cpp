#include<stdio.h>
#include<cmath>
#define PI 3.14159265358979323846

//Lax-Wendroff格式
double LaxWendroffScheme(double uj, double uj_1, double uj1, double dt, double dx) {
    return uj - dt / dx * (uj1 - uj_1) / 2 + dt * dt / dx / dx * (uj1 - 2 * uj + uj_1) / 2;
}

//Lax-Wendroff格式的L2误差
double LaxWendroffScheme_error(double dt, int n) {
    double dx;//dx空间步长
    dx = double(3) / double(n);
    double u1[n + 1], u2[n + 1];//创建两个u数组用于交替迭代

    //初始值
    for (int i = 0; i <= n; i++) {
        u1[i] = sin(2 * PI * i * dx);
    }

    int j;
    //计算1个时间单位， 自动计算循环次数
    for (j = 1; j <= round(double(1) / dt); j++) {
        //奇数次
        if (j % 2 != 0) {
            for (int i = 1; i < n; i++) {
                u2[i] = LaxWendroffScheme(u1[i], u1[i - 1], u1[i + 1], dt, dx);//用Lax-Wendroff格式迭代
            }
            //补充边界点值
            u2[0] = LaxWendroffScheme(u1[0], u1[n - 1], u1[1], dt, dx);
            u2[n] = u2[0];
        }
        //偶数次
        else {
            for (int i = 1; i < n; i++) {
                u1[i] = LaxWendroffScheme(u2[i], u2[i - 1], u2[i + 1], dt, dx);//用Lax-Wendroff格式迭代
            }
            //补充边界点值
            u1[0] = LaxWendroffScheme(u2[0], u2[n - 1], u2[1], dt, dx);
            u1[n] = u1[0];
        }
    }
    j--;

    double sum = 0;
    if (j % 2 != 0) {//如果最后存储在u2
        for (int i = 0; i <= n; i++) {
            sum += (u2[i] - sin(2 * PI * (i * dx - 1))) * (u2[i] - sin(2 * PI * (i * dx - 1)));
        }
    }
    else {//如果最后存储在u1
        for (int i = 0; i <= n; i++) {
            sum += (u1[i] - sin(2 * PI * (i * dx - 1))) * (u1[i] - sin(2 * PI * (i * dx - 1)));
        }
    }

    return sqrt(sum / n);
}

int main() {
    double dt;//dt初始时间步长
    int n;//初始空间网格等分数
    printf("请输入初始时间步长：");
    scanf("%lf", &dt);
    printf("请输入初始空间网格等分数：");
    scanf("%d", &n);
    double dx = double(3) / double(n);//初始网格间距
    double c = dt / dx;
    printf("初始空间步长：%lf, c:%lf", dx, c);

    double p;//精度
    while (n < 50000) {
        p = log2(LaxWendroffScheme_error(dt, n) / LaxWendroffScheme_error(dt / 2, n * 2));
        printf("dt=%lf，n=%d精度：%.4f\n", dt, n, p);
        dt = dt / 2;
        n = n * 2;
    }

    return 0;
}

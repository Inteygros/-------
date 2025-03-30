#include<stdio.h>
#include<cmath>
#define PI 3.14159265358979323846

//upwind Warming-Beam格式
double UpwindWarmingBeamScheme(double uj, double uj_1, double uj_2, double dt, double dx) {
    return uj - dt / dx * (uj - uj_1) - dt / dx / 2 * (1 - dt / dx) * (uj - 2 * uj_1 + uj_2);
}

//upwind Warming-Beam格式的L2误差
double UpwindWarmingBeamScheme_error(double dt, int n) {
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
            for (int i = 2; i <= n; i++) {
                u2[i] = UpwindWarmingBeamScheme(u1[i], u1[i - 1], u1[i - 2], dt, dx);//用upwind Warming-Beam格式迭代
            }
            //补充边界点值
            u2[1] = UpwindWarmingBeamScheme(u1[1], u1[0], u1[n - 1], dt, dx);
            u2[0] = u2[n];
        }
        //偶数次
        else {
            for (int i = 2; i <= n; i++) {
                u1[i] = UpwindWarmingBeamScheme(u2[i], u2[i - 1], u2[i - 2], dt, dx);//用upwind Warming-Beam格式迭代
            }
            //补充边界点值
            u1[1] = UpwindWarmingBeamScheme(u2[1], u2[0], u2[n - 1], dt, dx);
            u1[0] = u1[n];
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
        p = log2(UpwindWarmingBeamScheme_error(dt, n) / UpwindWarmingBeamScheme_error(dt / 2, n * 2));
        printf("dt=%lf，n=%d精度：%.4f\n", dt, n, p);
        dt = dt / 2;
        n = n * 2;
    }

    return 0;
}

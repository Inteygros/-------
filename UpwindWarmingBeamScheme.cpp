#include<stdio.h>
#include<cmath>
#define PI 3.14159265358979323846

//upwind Warming-Beam格式
double UpwindWarmingBeamScheme(double uj, double uj_1, double uj_2, double dt, double dx) {
    return uj - dt / dx * (uj - uj_1) - dt / dx / 2 * (1 - dt / dx) * (uj - 2 * uj_1 + uj_2);
}

int main() {
    double dt, dx, c;//dt时间步长，dx空间步长
    int n;//空间网格等分数
    printf("请输入时间步长：");
    scanf("%lf", &dt);
    printf("请输入空间网格等分数：");
    scanf("%d", &n);
    dx = double(3) / double(n);
    c = dt / dx;
    printf("空间步长：%lf,c:%lf", dx, c);
    double u1[n + 1], u2[n + 1];//创建两个u数组用于交替迭代

    // 定义文件名
    char filename[256];

    // 将dt和c的值写入文件名
    snprintf(filename, sizeof(filename), "UWBS_output_dt=%.4f_dx=%.4f.csv", dt, dx);

    FILE* output = fopen(filename, "w");//输出文件

    //初始值
    for (int i = 0; i <= n; i++) {
        u1[i] = sin(2 * PI * i * dx);
    }

    //计算10个时间单位， 自动计算循环次数
    for (int j = 1; j <= round(double(10) / dt); j++) {
        //奇数次
        if (j % 2 != 0) {
            for (int i = 2; i <= n; i++) {
                u2[i] = UpwindWarmingBeamScheme(u1[i], u1[i - 1], u1[i - 2], dt, dx);//用upwind Warming-Beam格式迭代
            }
            //补充边界点值
            u2[1] = UpwindWarmingBeamScheme(u1[1], u1[0], u1[n - 1], dt, dx);
            u2[0] = u2[n];

            //记录输出，记录误差
            for (int i = 0; i < n; i++) {
                fprintf(output, "%lf ", u2[i]);
            }
            fprintf(output, "%lf\n", u2[n]);
        }
        //偶数次
        else {
            for (int i = 2; i <= n; i++) {
                u1[i] = UpwindWarmingBeamScheme(u2[i], u2[i - 1], u2[i - 2], dt, dx);//用upwind Warming-Beam格式迭代
            }
            //补充边界点值
            u1[1] = UpwindWarmingBeamScheme(u2[1], u2[0], u2[n - 1], dt, dx);
            u1[0] = u1[n];

            //记录输出，记录误差
            for (int i = 0; i < n; i++) {
                fprintf(output, "%lf ", u1[i]);
            }
            fprintf(output, "%lf\n", u1[n]);
        }
    }

    fclose(output);
    printf("\n数据已保存到UWBS_output_dt=%.4f_dx=%.4f.csv\n", dt, dx);

    return 0;
}

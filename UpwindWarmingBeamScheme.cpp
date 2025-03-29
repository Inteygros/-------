#include<stdio.h>
#include<cmath>
#define PI 3.14159265358979323846

//一阶迎风格式
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
    double u1[n + 1], u2[n + 1], e[n + 1];//创建两个u数组用于交替迭代
    FILE* error = fopen("UWBS_error.csv", "w");//误差文件
    FILE* output = fopen("UWBS_output.csv", "w");//输出文件

    for (int i = 0; i <= n; i++) {
        u1[i] = sin(2 * PI * i * dx);
    }

    //计算100个时间步长，可修改
    for (int j = 1; j <= 100; j++) {
        //奇数次
        if (j % 2 != 0) {
            for (int i = 2; i <= n; i++) {
                u2[i] = UpwindWarmingBeamScheme(u1[i], u1[i - 1], u1[i - 2], dt, dx);//用upwind Warming-Beam格式迭代
                e[i] = (u2[i] - u1[i]) / dt + (3 * u1[i] - 4 * u1[i - 1] + u1[i - 2]) / 2 / dx - dt / 2 / dx / dx * (u1[i] - 2 * u1[i - 1] + u1[i - 2]);
            }
            //补充边界点值
            u2[1] = UpwindWarmingBeamScheme(u1[1], u1[0], u1[n - 1], dt, dx);
            e[1] = (u2[1] - u1[1]) / dt + (3 * u1[1] - 4 * u1[0] + u1[n - 1]) / 2 / dx - dt / 2 / dx / dx * (u1[1] - 2 * u1[0] + u1[n - 1]);
            u2[0] = u2[n];
            e[0] = e[n];

            //记录输出，记录误差
            for (int i = 0; i <= n; i++) {
                fprintf(output, "%.10e ", u2[i]);
                fprintf(error, "%.10e ", e[i]);
            }
            fprintf(output, "%.10e\n", u2[n]);
            fprintf(error, "%.10e\n", e[n]);
        }
        //偶数次
        else {
            for (int i = 2; i <= n; i++) {
                u1[i] = UpwindWarmingBeamScheme(u2[i], u2[i - 1], u2[i - 2], dt, dx);//用upwind Warming-Beam格式迭代
                e[i] = (u1[i] - u2[i]) / dt + (3 * u2[i] - 4 * u2[i - 1] + u2[i - 2]) / 2 / dx - dt / 2 / dx / dx * (u2[i] - 2 * u2[i - 1] + u2[i - 2]);
            }
            //补充边界点值
            u1[1] = UpwindWarmingBeamScheme(u2[1], u2[0], u2[n - 1], dt, dx);
            e[1] = (u1[1] - u2[1]) / dt + (3 * u2[1] - 4 * u2[0] + u2[n - 1]) / 2 / dx - dt / 2 / dx / dx * (u2[1] - 2 * u2[0] + u2[n - 1]);
            u1[0] = u1[n];
            e[0] = e[n];

            //记录输出，记录误差
            for (int i = 0; i <= n; i++) {
                fprintf(output, "%.10e ", u1[i]);
                fprintf(error, "%.10e ", e[i]);
            }
            fprintf(output, "%.10e\n", u1[n]);
            fprintf(error, "%.10e\n", e[n]);
        }
    }

    fclose(output);
    fclose(error);
    printf("\n数据已保存到UWBS_output.csv和UWBS_error.csv\n");
    return 0;
}

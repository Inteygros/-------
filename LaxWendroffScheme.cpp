#include<stdio.h>
#include<cmath>
#define PI 3.14159265358979323846

//Lax-Wendroff格式
double LaxWendroffScheme(double uj, double uj_1, double uj1, double dt, double dx) {
    return uj - dt / dx * (uj1 - uj_1) / 2 + dt * dt / dx / dx * (uj1 - 2 * uj + uj_1) / 2;
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
    snprintf(filename, sizeof(filename), "LWS_output_dt=%.4f_dx=%.4f.csv", dt, dx);

    FILE* output = fopen(filename, "w");//输出文件

    //初始值
    for (int i = 0; i <= n; i++) {
        u1[i] = sin(2 * PI * i * dx);
    }

    //计算10个时间单位， 自动计算循环次数
    for (int j = 1; j <= round(double(10) / dt); j++) {
        //奇数次
        if (j % 2 != 0) {
            for (int i = 1; i < n; i++) {
                u2[i] = LaxWendroffScheme(u1[i], u1[i - 1], u1[i + 1], dt, dx);//用Lax-Wendroff格式迭代
            }
            //补充边界点值
            u2[0] = LaxWendroffScheme(u1[0], u1[n - 1], u1[1], dt, dx);
            u2[n] = u2[0];

            //记录输出
            for (int i = 0; i < n; i++) {
                fprintf(output, "%lf ", u2[i]);
            }
            fprintf(output, "%lf\n", u2[n]);
        }
        //偶数次
        else {
            for (int i = 1; i < n; i++) {
                u1[i] = LaxWendroffScheme(u2[i], u2[i - 1], u2[i + 1], dt, dx);//用Lax-Wendroff格式迭代
            }
            //补充边界点值
            u1[0] = LaxWendroffScheme(u2[0], u2[n - 1], u2[1], dt, dx);
            u1[n] = u1[0];

            //记录输出，记录误差
            for (int i = 0; i < n; i++) {
                fprintf(output, "%lf ", u1[i]);
            }
            fprintf(output, "%lf\n", u1[n]);
        }
    }

    fclose(output);
    printf("\n数据已保存到LWS_output_dt=%.4f_dx=%.4f.csv\n", dt, dx);

    return 0;
}

#include<stdio.h>
#include<cmath>
#define PI 3.14159265358979323846

//一阶迎风格式
double FirstOrderUpwindScheme(double uj, double uj_1, double dt, double dx) {
    return uj - dt / dx * (uj - uj_1);
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
    double u[n + 1];

    // 定义文件名
    char filename[256];

    // 将dt和c的值写入文件名
    snprintf(filename, sizeof(filename), "FOUS_output_dt=%.4f_dx=%.4f.csv", dt, dx);

    FILE* output = fopen(filename, "w");//输出文件

    for (int i = 0; i <= n; i++) {
        u[i] = sin(2 * PI * i * dx);
    }

    //计算10个时间单位， 自动计算循环次数
    for (int j = 1; j <= round(double(10) / dt); j++) {
        for (int i = n; i > 0; i--) {
            double ui = u[i];
            u[i] = FirstOrderUpwindScheme(u[i], u[i - 1], dt, dx);
        }
        u[0] = u[n];

        //记录输出
        for (int i = 0; i < n; i++) {
            fprintf(output, "%.10e ", u[i]);
        }
        fprintf(output, "%.10e\n", u[n]);
    }

    fclose(output);
    printf("\n数据已保存到FOUS_output_dt=%.4f_dx=%.4f.csv\n", dt, dx);

    return 0;
}

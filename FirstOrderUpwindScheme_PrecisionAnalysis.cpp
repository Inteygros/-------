#include<stdio.h>
#include<cmath>
#define PI 3.14159265358979323846

//一阶迎风格式的L2误差
double FirstOrderUpwindScheme_error(double dt, int n) {
    double dx;//dt时间步长，dx空间步长
    dx = double(3) / double(n);
    double u[n + 1];

    for (int i = 0; i <= n; i++) {
        u[i] = sin(2 * PI * i * dx);
    }

    //计算1个时间单位， 自动计算循环次数
    for (int j = 1; j <= round(double(1) / dt); j++) {
        for (int i = n; i > 0; i--) {
            u[i] = u[i] - dt / dx * (u[i] - u[i - 1]);
        }
        u[0] = u[n];
    }

    //计算误差的平方和
    double sum = 0;
    for (int i = 0; i <= n; i++) {
        sum += (u[i] - sin(2 * PI * (i * dx - 10))) * (u[i] - sin(2 * PI * (i * dx - 10)));
    }

    return sqrt(sum / n);
}

int main() {
    double dt;//dt初始时间步长
    int n;//初始空间网格等分数
    printf("请输入时间步长：");
    scanf("%lf", &dt);
    printf("请输入空间网格等分数：");
    scanf("%d", &n);
    double dx = double(3) / double(n);//初始网格间距
    double c = dt / dx;
    printf("初始空间步长：%lf, c:%lf", dx, c);

    double p;//精度
    while (n < 50000) {
        p = log2(FirstOrderUpwindScheme_error(dt, n) / FirstOrderUpwindScheme_error(dt / 2, n * 2));
        printf("dt=%lf，n=%d精度：%.4f\n", dt, n, p);
        dt = dt / 2;
        n = n * 2;
    }

    return 0;
}
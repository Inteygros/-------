#include<stdio.h>
#include<cmath>
#define PI 3.14159265358979323846

//一阶迎风格式
double FirstOrderUpwindScheme(double uj, double uj_1, double a, double dt, double dx){
    return uj-a*dt/dx*(uj-uj_1);
}

int main(){
    double dt, dx, c;//dt时间步长，dx空间步长
    int n;//空间网格等分数
    printf("请输入时间步长：");
    scanf("%lf",&dt);
    printf("请输入空间网格等分数：");
    scanf("%d",&n);
    dx = double(3)/double(n);
    c = dt/dx;
    printf("空间步长：%lf,c:%lf",dx,c);
    double u[n+1],e[n+1];
    FILE *error = fopen("error.csv", "w");//误差文件
    FILE *output = fopen("output.csv", "w");//输出文件

    for(int i = 0; i <= n; i++){
        u[i]=sin(2 * PI * i * dx);
    }

    //计算100个时间步长，可修改
    for(int j = 1; j <= 100; j++){
        for(int i = n; i > 0; i--){
            double ui = u[i];
            u[i] = FirstOrderUpwindScheme(u[i], u[i-1], 1, dt, dx);
            e[i] = (u[i] - ui) / dt + (ui - u[i-1]) / dx;
        }
        u[0] = u[n];
        e[0] = e[n];

        //记录输出，记录误差
        for(int i = 0; i <= n; i++){
        }
        for(int i = 0; i < n; i++){
            fprintf(output,"%.10e ",u[i]);
            fprintf(error,"%.10e ",e[i]);
        }
        fprintf(output,"%.10e\n",u[n]);
        fprintf(error,"%.10e\n",e[n]);
    }

    fclose(output);
    fclose(error);
    printf("\n数据已保存到output.csv和error.csv\n");

    return 0;
}

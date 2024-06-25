#include <stdio.h>
#include <string.h>

int main() {
    FILE *src_file, *dst_file; // 创建指针变量，以存放文件
    char line[1024];
    int count = 0;

    src_file = fopen("result.txt", "r"); // 打开结果文件
    if (src_file == NULL) { // 如果结果文件不存在，执行以下代码
        printf("the file is not exist");
        return 1;
    }

    dst_file = fopen("destination.txt", "w"); // 写入文件
    if (dst_file == NULL) { // 如果写入文件不存在，执行以下代码
        printf("the file is not exist");
        fclose(src_file);
        return 1;
    }

    while (fgets(line, sizeof(line), src_file)) {
        if (strstr(line, "Forward:	Score:") != NULL) { // 在一个字符串中搜索另一个字符串
            count = 9;
            while (count >= 0 && fgets(line, sizeof(line), src_file)) {
                fputs(line, dst_file);
                count--;
            }
            // continue;
        }
    }

    fclose(src_file);
    fclose(dst_file);

    return 0;
}

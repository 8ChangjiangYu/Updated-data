#include <stdio.h>
#include <string.h>

int main() {
    FILE *src_file, *dst_file; // 
    char line[1024];
    int count = 0;

    src_file = fopen("result.txt", "r"); // 
    if (src_file == NULL) { // 
        printf("the file is not exist");
        return 1;
    }

    dst_file = fopen("destination.txt", "w"); // 
    if (dst_file == NULL) { // 
        printf("the file is not exist");
        fclose(src_file);
        return 1;
    }

    while (fgets(line, sizeof(line), src_file)) {
        if (strstr(line, "Forward:	Score:") != NULL) { // 
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

#include <stdio.h>

int main()
{
  int x = 1;

  char *y = (char*)&x;

  // printf("If 1: little endian; if 0: big endian\n");
  printf("%c\n",*y+48);

  return  0;
}

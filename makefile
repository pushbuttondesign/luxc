all:
	gcc -ggdb -Wall -Wextra -O0 -std=c11 \
	-o ./example.out \
	./example.c \
	../SPDlibrary/SPDlib.c \
	-I /root/code/C/librarys/include \
	-lm

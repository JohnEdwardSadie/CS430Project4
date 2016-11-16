all: main.c
	gcc -o raytrace main.c

clean:
	rm -rf raycast *~


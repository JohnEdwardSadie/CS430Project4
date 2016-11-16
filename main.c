#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


//John E. Sadie
//CS430 Computer Graphics
//Project 4 - Recursive Raytracing
//11.15.16

/*
Description:
In the previous project you will wrote code to raycast and shade mathematical primitives based
on a scene input file into a pixel buffer. In this project you will add recursive raytracing to
provide reflection and refraction.
*/


//declaration list of variables
//From the JSON file
typedef struct{
        char type;
        double color[3];
        double diffusecolor[3];
        double specularColor[3];
        double position[3];
        double normal[3];
        double radius;
        double reflectivity;
        double refractivity;
        double refractionIndex;
        double width, height;
        unsigned char r,g,b;
} Scene;

typedef struct{
        char type;
        double color[3];
        double position[3];
        double direction[3];
        double theta;
        double radiala2;
        double radiala1;
        double radiala0;
        double angulara0;
} Light;

Scene r,g,b;
Scene camera;
Scene *scene;
Light *lightScene;
Scene *PixelBuffer;

int line = 1;
int cameraOne;
int incrementObject;
int lastIndex = 0;
int lastIndexLight = 0;
Pi = 3.14;
double black[3] = {0,0,0};

/*Function declaration for shading*/
double* Shade(int recurseValue, double* color, int closestObject, double Ro[], double Rd[], double closestT, double rIndex);
/*Math functions*/

//creating a square operator
static inline double Sqr(double n){
	return n*n;
}

static inline void normalize(double* v){
    //getting the length of the vector
    double length = sqrt(Sqr(v[0]) + Sqr(v[1]) + Sqr(v[2]));

    //normalizing our vector
    //turning it into a unit vector
    v[0] = v[0]/length;
    v[1] = v[1]/length;
    v[2] = v[2]/length;
}
/*The following equations were provided by Palmer in-class(hand-written)*/
//Sources: http://csis.pace.edu/~marchese/CG_Rev/Lect10New/cg_l10new.htm
static inline double VectorLength(double* v){
    return sqrt(Sqr(v[0])+Sqr(v[1])+Sqr(v[2]));

}
static inline void VectorScaling(double scaler, double* vector, double* ans){
    ans[0] = scaler * vector[0];
    ans[1] = scaler * vector[1];
    ans[2] = scaler * vector[2];
}

static inline void VectorSubtraction(double* a, double* b, double* ans){
    ans[0] = a[0] - b[0];
    ans[1] = a[1] - b[1];
    ans[2] = a[2] - b[2];
}

static inline void VectorAddition(double* a, double* b, double* ans){
    ans[0] = a[0] + b[0];
    ans[1] = a[1] + b[1];
    ans[2] = a[2] + b[2];
}

//clamps the color(restrict the value)
static inline double clamp(double n){
    if(n > 1){
        n = 1;
    }
    else if(n < 0){
        n = 0;
    }
    else{
        return n;
    }
}

static inline double DotProduct(double* a, double* b){
    return ((a[0]*b[0] + a[1]*b[1] + a[2]*b[2]));
}
static inline void VectorReflection(double* normal, double* R, double* reflection){
    double a[3] = {0,0,0};
    double nv = DotProduct(normal, R);

    VectorScaling(nv,normal, a);
    VectorScaling(2,a,a);
    VectorSubtraction(R,a,reflection);
}

//Source: http://csis.pace.edu/~marchese/CG_Rev/Lect10New/cg_l10new.htm
static inline double AngularAttenuation(Light a, double* NewRd, double Pi){
    double SpotlightToObject[3];
    SpotlightToObject[0] = (1 * NewRd[0]);
    SpotlightToObject[1] = (1 * NewRd[1]);
    SpotlightToObject[2] = (1 * NewRd[2]);
    double cos = (((a.direction[0]*SpotlightToObject[0]) +(a.direction[1]*SpotlightToObject[1]) + (a.direction[2]*SpotlightToObject[2])));
    return pow(cos, a.angulara0);
    free(SpotlightToObject);
}

//Source: http://csis.pace.edu/~marchese/CG_Rev/Lect10New/cg_l10new.htm
static inline double RadialAttenuation(Light a, double* NewRo){
    double* VL = malloc(sizeof(double)*3);
    VectorSubtraction(a.position, NewRo, VL);
    double den = (a.radiala2 * Sqr(VectorLength(VL))) + (a.radiala1 * VectorLength(VL)) + (a.radiala0);
    return 1/den;
    free(VL);
}
//Source: http://csis.pace.edu/~marchese/CG_Rev/Lect10New/cg_l10new.htm
static inline double IlluminateDiffuse(int index,double* diffusecolor, Light objectlight, double NL){
    if(NL <= 0) return 0;
    return diffusecolor[index]*objectlight.color[index]*NL;
}
/*
   Calculate normal vector using cross product:
n = PxQ
nx = PyQz - PzQy
ny = PzQx -PxQz
nz = PxQy - PyQx
Create a light vector
�       Select a point on the polygon surface and position for the light source
     A polygon vertex may be selected or a position inside the polygon.
*/
//Source: http://csis.pace.edu/~marchese/CG_Rev/Lect10New/cg_l10new.htm
static inline double IlluminateSpecular(int index, double* SpecularColor, Light objectlight, double NL, double VR, double ns){
    if(VR <= 0){
        return 0;
    }
    if(NL <= 0){
        return 0;
    }

    return SpecularColor[index]*objectlight.color[index]*pow(VR,ns);
}

//implements pseudocode provided by Palmer
double sphereIntersection(double* Ro, double* Rd, double* position, double radius){
    double a, b, c;
    normalize(Rd);

    //sphere intersection equation provided by Dr. Palmer
    a = Sqr(Rd[0])+Sqr(Rd[1])+Sqr(Rd[2]);
    b = 2*(Rd[0]*(Ro[0]-position[0]) + Rd[1]*(Ro[1]-position[1]) + Rd[2]*(Ro[2]-position[2]));
    c = (Sqr((Ro[0]-position[0])) + Sqr((Ro[1]-position[1])) + Sqr((Ro[2]-position[2])) - Sqr(radius));

    double t0 = ((-b - sqrt(Sqr(b) - 4.0*c*a))/(2.0*a));
    double t1 = ((-b + sqrt(Sqr(b) - 4.0*c*a))/(2.0*a));


    if(t0 > 0.0){
        return t0;
    }
    if(t1 > 0.0){
        return t1;
    }
    return -1;

    double desc = (Sqr(b) - 4*a*c);

    if(desc < 0.0){
        return 999999;
    }
}

//planeIntersection function
//Source: http://www.siggraph.org/education/materials/HyperGraph/raytrace/rayplane_intersection.htm
double planeIntersection(double* Ro, double* Rd, double* position, double* normal){
    normalize(normal);
    normalize(Rd);

    //The length from camera to plane
    double d = -(normal[0]*position[0] + normal[1]*position[1] + normal[2]*position[2]);

    //denominator
    double denominator = (normal[0]*Rd[0] + normal[1]*Rd[1] + normal[2]*Rd[2]);

    //plane intersection equation provided by Dr. Palmer
    //t = -(AX0 + BY0 + CZ0 + D) / (AXd + BYd + CZd)
    if(denominator == 0.0){
        return -1;
    }
    double t = -(normal[0]*Ro[0] + normal[1]*Ro[1] + normal[2]*Ro[2] + d)/(normal[0]*Rd[0] + normal[1]*Rd[1] + normal[2]*Rd[2]);

    return t;
}
//RayCast function
void rayCast(double N, double M){
    //instantiating variables
    int i = 0;
    int j = 0;
    int index;
    double Ro[3] = {0, 0, 0};
    double center[3] = {0, 0, 0};


    double w = camera.width;
    double h = camera.height;
    double PixelWidth = w/N;
    double PixelHeight = h/M;

    double p_z = 1;
    for(i=0; i<M; i+=1){
        for(j=0; j<N; j+=1){
            double p_y = center[1] - h/2.0 + PixelHeight*(i+0.5);
            double p_x = center[0] - w/2.0 + PixelWidth*(j+0.5);
            double p_z = 1; //z-coordinate view plane
            double Rd[3] = {p_x, p_y, p_z};
            //Normalization of Rd
            normalize(Rd);

            double closestT = 999999; //closest point to the camera
            int closestObject = -1;

            for(index=0; index<=lastIndex; index++){

                double t = 0;
                //Shoot function from the pseudocode provided by Palmer
                if(scene[index].type == 's'){
                        t = sphereIntersection(Ro, Rd, scene[index].position, scene[index].radius);
                }
                //Shoot function from the pseudocode provided by Palmer
                if(scene[index].type == 'p'){
                        t = planeIntersection(Ro, Rd, scene[index].position, scene[index].normal);

                }
                if(t > 0 && t < closestT){
                    closestT = t;
                    closestObject = index;
                }
            }
            double color[3] = {0,0,0};
            color[0] = 0;
            color[1] = 0;
            color[2] = 0;

            double actual[3] = {0,0,0};
            if(closestT > 0 && closestT != 999999){
                memcpy(actual, Shade(7, color, closestObject, Ro, Rd, closestT, 1), sizeof(double)*3);
            }
            int pos = (int)((M - i - 1)*N + j);
            PixelBuffer[pos].r = (char)(clamp(actual[0])*255);
            PixelBuffer[pos].g = (char)(clamp(actual[1])*255);
            PixelBuffer[pos].b = (char)(clamp(actual[2])*255);
        }
    }
}

//Shade function
double* Shade(int recurseValue, double* color, int closestObject, double Ro[], double Rd[], double closestT, double rIndex){
    int i;
    int j;
    //Setting the base case for recursion
    if(recurseValue == 0 || closestObject < 0 || closestT == 999999){
        return black;
    }
    double NewRo[3] = {0,0,0};
    NewRo[0] = (closestT * Rd[0]) + Ro[0];
    NewRo[1] = (closestT * Rd[1]) + Ro[1];
    NewRo[2] = (closestT * Rd[2]) + Ro[2];

    /* Ambient + Diffuse + Emission Light */
    for(i=0; i<=lastIndexLight; i++){
        double NewRd[3] = {0,0,0};
        VectorSubtraction(lightScene[i].position, NewRo, NewRd);
        normalize(NewRd);
        double DP = sqrt(Sqr(lightScene[i].position[0] - NewRo[0])+ Sqr(lightScene[i].position[1] - NewRo[1])+ Sqr(lightScene[i].position[2] - NewRo[2]));
        double closestShadow = 999999;
        int closestSIndex = -1;
        //Checks what casts a shadow
        for(j=0; j<=lastIndex; j++){
            double closestS = 0;
            if(j == closestObject){
                continue;
            }
            if(scene[j].type == 's'){
                closestS = sphereIntersection(NewRo, NewRd, scene[j].position, scene[j].radius);
            }
            if(scene[j].type == 'p'){
                closestS = planeIntersection(NewRo, NewRd, scene[j].position, scene[j].normal);
            }
            if(closestT > DP){
                continue;
            }
            if(closestS > 0.0 && closestS < closestShadow){
                closestShadow = closestS;
                closestSIndex = j;
            }
        }
        if(closestSIndex < 0){
            double closestN[3] = {0,0,0};
            if(scene[closestObject].type == 's') {
                    VectorSubtraction(NewRo, scene[closestObject].position, closestN);
                    }
            if(scene[closestObject].type == 'p'){
                closestN[0] = scene[closestObject].normal[0];
                closestN[1] = scene[closestObject].normal[1];
                closestN[2] = scene[closestObject].normal[2];
            }
            //light position from intersection
            normalize(closestN); //grab real position
            double *DL = NewRd;
            double reflect[3] = {0,0,0};
            VectorReflection(closestN, DL, reflect);
            normalize(reflect);
            //direction of the object seen by camera
            double cameraToObject[3] = {0,0,0};
            cameraToObject[0] = Rd[0];
            cameraToObject[1] = Rd[1];
            cameraToObject[2] = Rd[2];
            normalize(cameraToObject);

            //current diffuse and specular
            double* Diffuse = scene[closestObject].diffusecolor;
            double* Specular = scene[closestObject].specularColor;
            //Palmer suggested 20
            //For shininess
            double ns = 20;
            normalize(closestN);
            normalize(DL);
            double VR = DotProduct(cameraToObject, reflect);
            double NL = DotProduct(closestN, DL);

            //Values of reflection
              double recurseRo[3] = {0,0,0};
              memcpy(recurseRo,NewRo, sizeof(double)*3);

              double recurseRd[3] = {0,0,0};
              double c1;
              c1 = DotProduct(closestN, Rd);
              c1 = -1*c1;
              //Using our math functions instantiated in the beginning for scaling
              VectorScaling(c1, closestN, recurseRd);
              VectorScaling(2, recurseRd, recurseRd);
              VectorAddition(Rd, recurseRd, recurseRd);
              normalize(recurseRd);

              recurseRo[0] = NewRo[0] + recurseRd[0]*0.01;
              recurseRo[1] = NewRo[1] + recurseRd[1]*0.01;
              recurseRo[2] = NewRo[2] + recurseRd[2]*0.01;



              double reflection_color[3] = {0,0,0};
              double refraction_color[3] = {0,0,0};
              if(scene[closestObject].reflectivity != 0){
                  Shade(recurseValue-1, reflection_color, RI, recurseRo, recurseRd, closestRT,1);
              }
              /* Implementing the color of reflectivity */

              double e = (1 - scene[closestObject].reflectivity - scene[closestObject].refractivity);
              //0
              color[0] +=  e * (RadialAttenuation(lightScene[i], NewRo) * AngularAttenuation(lightScene[i], NewRd, Pi)) *
                  ((IlluminateDiffuse(0, Diffuse, lightScene[i], NL) + IlluminateSpecular(0, Specular, lightScene[i], NL, VR, ns))) +
                  (((scene[closestObject].reflectivity) * reflection_color[0]));
              //1
              color[1] += e * (RadialAttenuation(lightScene[i], NewRo) * AngularAttenuation(lightScene[i], NewRd, Pi)) *
                  ((IlluminateDiffuse(1, Diffuse, lightScene[i], NL) + IlluminateSpecular(1,Specular,lightScene[i], NL, VR, ns))) +
                  (((scene[closestObject].reflectivity) * reflection_color[1]));
              //2
              color[2] += e * (RadialAttenuation(lightScene[i], NewRo) * AngularAttenuation(lightScene[i],NewRd,Pi)) *
                  ((IlluminateDiffuse(2, Diffuse, lightScene[i], NL) + IlluminateSpecular(2, Specular, lightScene[i], NL, VR, ns))) +
                  (((scene[closestObject].reflectivity) * reflection_color[2]));

        }
    }
    return color;
}



// next_c() wraps the getc() function and provides error checking and line
// nvber maintenance
int next_c(FILE* json) {
  int c = fgetc(json);
#ifdef DEBUG
  printf("next_c: '%c'\n", c);
#endif
  if (c == '\n') {
    line += 1;
  }
  if (c == EOF) {
    fprintf(stderr, "Error: Unexpected end of file on line nvber %d.\n", line);
    exit(1);
  }
  return c;
}


// expect_c() checks that the next character is d.  If it is not it emits
// an error.
void expect_c(FILE* json, int d) {
  int c = next_c(json);
  if (c == d) return;
  fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, line);
  exit(1);
}


// skip_ws() skips white space in the file.
void skip_ws(FILE* json) {
  int c = next_c(json);
  while (isspace(c)) {
    c = next_c(json);
  }
  ungetc(c, json);
}


// next_string() gets the next string from the file handle and emits an error
// if a string can not be obtained.
char* next_string(FILE* json) {
  char buffer[129];
  int c = next_c(json);
  if (c != '"') {
    fprintf(stderr, "Error: Expected string on line %d.\n", line);
    exit(1);
  }
  c = next_c(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: Strings longer than 128 characters in lengthgth are not supported.\n");
      exit(1);
    }
    if (c == '\\') {
      fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
      exit(1);
    }
    if (c < 32 || c > 126) {
      fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
      exit(1);
    }
    buffer[i] = c;
    i += 1;
    c = next_c(json);
  }
  buffer[i] = 0;
  return strdup(buffer);
}

double next_number(FILE* json) {
  double value;
  fscanf(json, "%lf", &value);
  //printf("Value is: %lf\n", value);
  // Error check this..
  return value;
}

double* next_vector(FILE* json) {
  double* v = malloc(3*sizeof(double));
  expect_c(json, '[');
  skip_ws(json);
  v[0] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[1] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[2] = next_number(json);
  skip_ws(json);
  expect_c(json, ']');
  return v;
}

//read_scene function
//Read's in our json file
//populates our object array
//saves data for future usage
void read_scene(char* filengthame) {
    int index = -1;
    int c;
    char Object;

    FILE* json = fopen(filengthame, "r");

  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", filengthame);
    exit(1);
  }

  skip_ws(json);

  // Find the beginning of the list
  expect_c(json, '[');

  skip_ws(json);

  // Find the incrementObject

  while (1) {
    c = fgetc(json);
    if (c == ']') {
       printf("%d", line);
      fprintf(stderr, "Error: This is the worst scene file EVER.\n");
      fclose(json);
      return;
    }
    if (c == '{') {
      skip_ws(json);

      // Parse the object
    char* key = next_string(json);
    if (strcmp(key, "type") != 0) {

        fprintf(stderr, "Error: Expected \"type\" key on line nvber %d.\n", line);
        exit(1);
      }
    skip_ws(json);
    expect_c(json, ':');
    skip_ws(json);
    char* value = next_string(json);
    if (strcmp(value, "camera") == 0) {
        cameraOne += 1;
        Object = 'c';
      }
      //What object are we looking at
      //String compare of value and sphere
    else if (strcmp(value, "sphere") == 0) {
        incrementObject += 1;
        Object = 's';
        scene[lastIndex].type = 's';

      }
      //String compare of value and plane
    else if (strcmp(value, "plane") == 0) {
        incrementObject += 1;
        Object = 'p';
        scene[lastIndex].type = 'p';
      }
      //String compare of value and light
      else if (strcmp(value, "light") == 0) {
        incrementObject += 1;
        Object = 'l';
        lightScene[lastIndexLight].type = 'l';
      }

    else {
        fprintf(stderr, "Error: Unknown type, \"%s\", on line nvber %d.\n", value, line);
        exit(1);
      }
    skip_ws(json);
    int incrementCamera = 0;
    while (1) {
	// , }

        c = next_c(json);
        if (c == '}') {
            // stop parsing this object
            break;
        }
        else if (c == ',') {
            // read another field
            skip_ws(json);
            char* key = next_string(json);
            skip_ws(json);
            expect_c(json, ':');
            skip_ws(json);
            if ((strcmp(key, "width") == 0)         ||
                (strcmp(key, "height") == 0)        ||
                (strcmp(key, "radius") == 0)        ||
                (strcmp(key, "radial-a2") == 0)     ||
                (strcmp(key, "radial-a1") == 0)     ||
                (strcmp(key, "radial-a0") == 0)     ||
                (strcmp(key, "angular-a0") == 0)    ||
                (strcmp(key, "theta") == 0)         ||
                (strcmp(key, "ior") == 0)           ||
                (strcmp(key, "reflectivity") == 0)  ||
                (strcmp(key, "refractivity") == 0))
                {
                double value = next_number(json);
                //*populating object array with our json contents*//
                if((strcmp(key, "width") == 0)){
                    camera.width = value;
                    incrementCamera+=1;
                }
                else if((strcmp(key, "height") == 0)){
                    camera.height = value;
                    incrementCamera +=1;
                }
                else if((strcmp(key, "radius") == 0)){
                    scene[lastIndex].radius = value;
                }
                else if((strcmp(key, "radial-a2") == 0)){
                    lightScene[lastIndexLight].radiala2 = value;
                }
                else if((strcmp(key, "radial-a1") == 0)){
                    lightScene[lastIndexLight].radiala1 = value;
                }
                else if((strcmp(key, "radial-a0") == 0)){
                    lightScene[lastIndexLight].radiala0 = value;
                }
                else if((strcmp(key, "angular-a0") == 0)){
                    lightScene[lastIndexLight].angulara0 = value;
                }
                else if((strcmp(key, "theta") == 0)){
                    lightScene[lastIndexLight].theta = value;
                }
                else if((strcmp(key, "reflectivity") == 0)){
                    scene[lastIndex].reflectivity = value;
                }
                else if((strcmp(key, "refractivity") == 0)){
                    scene[lastIndex].refractivity = value;
                }
                else if((strcmp(key, "ior") == 0)){
                    scene[lastIndex].refractionIndex = value;
                }
            }
            else if ((strcmp(key, "color") == 0) ||
                     (strcmp(key, "position") == 0) ||
                     (strcmp(key, "normal") == 0) ||
                     (strcmp(key, "direction") == 0) ||
                     (strcmp(key, "diffuse_color") == 0) ||
                     (strcmp(key, "specular_color") == 0))
                        {
                double* value = next_vector(json);
                //assigning value of color
                if((strcmp(key, "color") == 0)){
                        //set color value of scene
                    if(Object != 'l'){
                        scene[lastIndex].color[0] = value[0];
                        scene[lastIndex].color[1] = value[1];
                        scene[lastIndex].color[2] = value[2];

                    }
                    //set color value of lightScene
                    else if(Object == 'l'){
                        lightScene[lastIndexLight].color[0] = value[0];
                        lightScene[lastIndexLight].color[1] = value[1];
                        lightScene[lastIndexLight].color[2] = value[2];
                    }
                }
                 //assigning value of diffusecolor
                else if((strcmp(key, "diffuse_color") == 0)){
                    scene[lastIndex].diffusecolor[0] = value[0];
                    scene[lastIndex].diffusecolor[1] = value[1];
                    scene[lastIndex].diffusecolor[2] = value[2];

                }
                //assigning value of specularColor
               else if((strcmp(key, "specular_color") == 0)){
                    scene[lastIndex].specularColor[0] = value[0];
                    scene[lastIndex].specularColor[1] = value[1];
                    scene[lastIndex].specularColor[2] = value[2];

                }
                //assigning value of position
               else if((strcmp(key, "position") == 0)){
                    //setting value of positon in scene
                    if(Object != 'l'){
                        scene[lastIndex].position[0] = value[0];
                        scene[lastIndex].position[1] = value[1];
                        scene[lastIndex].position[2] = value[2];

                    }
                    //setting value of of position in lightScene
                   else if(Object == 'l'){
                        lightScene[lastIndexLight].position[0] = value[0];
                        lightScene[lastIndexLight].position[1] = value[1];
                        lightScene[lastIndexLight].position[2] = value[2];

                    }
                }
                //assigning value of normal
               else if((strcmp(key, "normal") == 0)){
                    scene[lastIndex].normal[0] = value[0];
                    scene[lastIndex].normal[1] = value[1];
                    scene[lastIndex].normal[2] = value[2];

                }
                //assigning value of direction
               else if((strcmp(key, "direction") == 0)){
                    lightScene[lastIndexLight].direction[0] = value[0];
                    lightScene[lastIndexLight].direction[1] = value[1];
                    lightScene[lastIndexLight].direction[2] = value[2];
                }
            }
            else {
                fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n",
                key, line);

            }
            skip_ws(json);
        }
        else {
            fprintf(stderr, "Error: Unexpected value on line %d\n", line);
            exit(1);
        }
      }
        skip_ws(json);
    c = next_c(json);
    if (c == ',') {
        //Error checking
        skip_ws(json);
        if(Object != 'c'){
            if(Object != 'l'){
                lastIndex++;
            }
            if(Object == 'l'){
                lastIndexLight++;
            }
        }
        if(Object == 'c'){
            if(incrementCamera != 2){
                exit(1);
            }
        incrementCamera = 0;
        }
    }
    else if (c == ']') {
        //Iterated through all objects
        fclose(json);
        return;
    }
    else {
        fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
        exit(1);
      }
    }
  }
}

/*This was to testing to see if I've properly acquired
the contents within the JSON by populating my array of incrementObject.
I then printed it out to the terminal.
*/

/*
void printScene(){
    int index = 0;
    while(scene.object[index].color != NULL){
    printf("object: %d\n", index);
    printf("type %s\n", scene.object[index].type);
    printf("color: %f %f %f\n", scene.object[index].color[0],scene.object[index].color[1],scene.object[index].color[2]);
    printf("position: %f %f %f\n", scene.object[index].position[0],scene.object[index].position[1],scene.object[index].position[2]);
    if(scene.object[index].normal != NULL){
    printf("normal: %f %f %f\n", scene.object[index].normal[0],scene.object[index].normal[1],scene.object[index].normal[2]);
    }
    else{
    printf("radius: %f\n", scene.object[index].radius);
    }
    index++;
    }
}*/


//Customized write function from project 1
int write(int w, int h, FILE* outputFile){
    int i;

    FILE *fp;

    char BufferSize[2] = {'P', '6'};
    int height = h;
    int width = w;

    fp = fopen(outputFile, "wb");
    fwrite(BufferSize, sizeof(BufferSize), sizeof(BufferSize)-1, fp);

    fprintf(fp,"\n%d %d", width, height);
    fprintf(fp,"\n%d", 255);
    fprintf(fp,"\n");
            for (i=0; i < width*height; i++){
                fwrite(&PixelBuffer[i].r, 1, 1, fp);
                fwrite(&PixelBuffer[i].g, 1, 1, fp);
                fwrite(&PixelBuffer[i].b, 1, 1, fp);
            }
    return 0;
}


//Main function
int main(int argc, char** argv) {

    //memory allocation
    scene = malloc(sizeof(Scene)*128);
    lightScene = malloc(sizeof(Light)*128);

    //ascii to integer.
    double N = (double)atoi(argv[1]);
    double M = (double)atoi(argv[2]);
    //referencing the type for memory allocation
    PixelBuffer = (Scene*)malloc(sizeof(Scene)* N * M);
    //treating PixelBuffer as a series of bytes
    memset(PixelBuffer, 0, 3*N * M);

    read_scene(argv[3]);
    rayCast(N, M);

    //Utilizing our command line arguments
    write(N, M, argv[4]);

    return 0;
}

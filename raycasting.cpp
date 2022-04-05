#include <iostream>
#include <fstream>
#include <string>
#include <math.h>	//For sqrt() and acos()
#include <stdlib.h>

//PPM (Portable PixMap) format http://netpbm.sourceforge.net/doc/ppm.html
//or more palatably, https://en.wikipedia.org/wiki/Netpbm
//Has support for 0-255 colours in each RGB channel; 24 bits per pixel


//************************************** Type Declaration **************************************//

// Defined struct for point/vector in 3D space 
struct vec3 {
  float x;
  float y;
  float z;
};

//Defined struct for representing the values of a sphere
struct SphereType {
  vec3 *c;		// center
  float r;	        // radius
  int m = -1;		// material color
};

//Defined struct for representing the values of a cylinder
struct CylinderType {
  vec3 *c;       	// center
  vec3 *dir;	// direction
  float r;	// radius
  float len;	// length
};

struct ColorType {
  float r;
  float g;
  float b;
};

//Define each pixel as a struct
struct Pixel{
  int r;	// red
  int g;	// green
  int b;	// blue
};

//Define each ray as a struct
struct RayType {
  vec3 *p;      // point
  vec3 *dir;	// direction
};

//Define mtlcolor
struct MTLColorType{
  float odr;
  float odg;
  float odb;
  float osr;
  float osg;
  float osb;
  float ka;
  float kd;
  float ks;
  float n;
};

//Define light source
struct LightType{
  vec3 *location;
  float w;
  vec3 *pix;
};

//Define Light Source attenuation
struct AttenLightType{
  LightType *lig;
  float c1;
  float c2;
  float c3;
};

struct DepthCueType{
  int dcr;
  int dcg;
  int dcb;
  int aMax;
  int aMin;
  int distMax;
  int distMin;
};



//************************************** Helper Functions **************************************//


//********** Vector Helper Functions **********//

// for vectors to simplify modifying a vec3 variable
vec3 * scaler(vec3 *a, float scale) {
  vec3 *temp = (vec3*) calloc(1,sizeof(vec3));
  temp->x = a->x * scale;
  temp->y = a->y * scale;
  temp->z = a->z * scale;
  return temp;
}

// for finding the length of vec3 vector
float lengthVec3(vec3 *a) {
  return sqrt((a->x * a->x) + (a->y * a->y) + (a->z * a->z));
}

// for finding the unit vector of vec3
vec3 * unitVec3(vec3 *a) {
  vec3 *temp = (vec3*) calloc(1,sizeof(vec3));
  float len = lengthVec3(a);
  temp->x = (int) a->x / len;
  temp->y = a->y / len;
  temp->z = a->z / len;
  return temp;
}

// for adding two vec3 (a + b)
vec3 * addVec3(vec3 *a, vec3 *b) {
  vec3 *temp = (vec3*) calloc(1,sizeof(vec3));
  temp->x = a->x + b->x;
  temp->y = a->y + b->y;
  temp->z = a->z + b->z;
  return temp;
}

// for subtracting two vec3 (a - b)
vec3 * subVec3(vec3 *a, vec3 *b) {
  vec3 *temp = (vec3*) calloc(1,sizeof(vec3));
  temp->x = a->x - b->x;
  temp->y = a->y - b->y;
  temp->z = a->z - b->z;
  return temp;
}

// for calculating dot product between two vectors
float dotVec3(vec3 *a, vec3 *b) {
  return (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
}

/*
// for calculating the angle between two vectors
float angleCosVec3(vec3 a, vec3 b) {
	return acos(dotVec3(a, b) / (lengthVec3(a) * lengthVec3(b)));
}

// orthogonal projection of vectors (project a onto b)
vec3 projVec3(vec3 a, vec3 b) {
	float temp = dotVec3(a, b) / dotVec3(b, b);
	return scaler(b, temp);
}

// orthogonal projection of vectors but calculating the length of the new vector
float projLengthVec3(vec3 a, vec3 b) {
	return dotVec3(a, b) / lengthVec3(b);
}
*/

// returns the cross product of two vectors (a X b != b X a)
vec3 * crossProdVec3(vec3 *a, vec3 *b) {
  vec3 *temp = (vec3*) calloc(1,sizeof(vec3));
  temp->x = a->y * b->z - a->z * b->y;
  temp->y = a->z * b->x - a->x * b->z;
  temp->z = a->x * b->y - a->y * b->x;
  return temp;
}

/*
float angleSinVec3(vec3 a, vec3 b) {
	return asin(crossProdVec3(a, b) / (lengthVec3(a) * lengthVec3(b)));
}
*/

// determines if two vectors are the same
bool equalVec3(vec3 *a, vec3 *b) {
  if (a->x == b->x && a->y == b->y && a->z == b->z)
    return true;
  else return false;
}




//********** Ray and Color Helper Functions **********//

//similar set-up to trace_ray
//calculate shadow-ray and check for intersections between light source in ray/sphere intersection
float shadow_intersection(vec3 *lightLoc, float light_w, vec3 *origin, SphereType sphere[], int sphereCount, int id) {
  //determine the direction from ray/sphere to light source
  vec3 *dir = subVec3(lightLoc, origin);
  dir = unitVec3(dir);

  float closest_t = 999999999.0;
  int sphere_id = -1;
  
  //use quadratic formula to calculate intersections between light and ray/sphere point
  for (int i = 0; i < sphereCount; i++) {
    float B = (dir->x * (origin->x - sphere[i].c->x));
    B += (dir->y * (origin->y - sphere[i].c->y));
    B += (dir->z * (origin->z - sphere[i].c->z));
    B *= 2;
    float C = (origin->x - sphere[i].c->x) * (origin->x - sphere[i].c->x);
    C += ((origin->y - sphere[i].c->y) * (origin->y - sphere[i].c->y));
    C += ((origin->z - sphere[i].c->z) * (origin->z - sphere[i].c->z));
    C -= (sphere[i].r * sphere[i].r);
    
    float discriminant = (B * B) - (4 * C);
    
    if (discriminant > 0) {
      float t_pos = (float)((-1 * B) + sqrt(discriminant)) / 2;
      float t_neg = (float)((-1 * B) - sqrt(discriminant)) / 2;
      
      //positional light source; take distance into consideration
      if (light_w == 1.0) {
	if (sphere_id == -1) {
	  if (t_pos > 0 && t_neg > t_pos && i != id) {
	    vec3 *temp = scaler(dir, t_pos);
	    temp = addVec3(temp, origin);
	    
	    vec3 *temp_l = subVec3(lightLoc, origin);
	    if (lengthVec3(temp) < lengthVec3(temp_l)) {
	      closest_t = t_pos;
	      sphere_id = i;
	    }
	  }
	  else if (t_neg > 0 && t_pos > t_neg && i != id) {
	    vec3 *temp = scaler(dir, t_neg);
	    temp = addVec3(temp, origin);
	    
	    vec3 *temp_l = subVec3(lightLoc, origin);
	    if (lengthVec3(temp) < lengthVec3(temp_l)) {
	      closest_t = t_neg;
	      sphere_id = i;
	    }
	  }
	}
	else {
	  if (t_pos > 0 && t_neg > t_pos && closest_t > t_pos && i != id) {
	    vec3 *temp = scaler(dir, t_pos);
	    temp = addVec3(temp, origin);
	    
	    vec3 *temp_l = subVec3(lightLoc, origin);
	    if (lengthVec3(temp) < lengthVec3(temp_l)) {
	      closest_t = t_pos;
	      sphere_id = i;
	    }
	  }
	  else if (t_neg > 0 && t_pos > t_neg && closest_t > t_neg && i != id) {
	    vec3 *temp = scaler(dir, t_neg);
	    temp = addVec3(temp, origin);
	    
	    vec3 *temp_l = subVec3(lightLoc, origin);
	    if (lengthVec3(temp) < lengthVec3(temp_l)) {
	      closest_t = t_neg;
	      sphere_id = i;
	    }
	  }
	}
      }
      
      //directional source; calculate intersections normally
      else {
	if (sphere_id == -1) {
	  if (t_pos > 0 && t_neg > t_pos && i != id) {
	    closest_t = t_pos;
	    sphere_id = i;
	  }
	  else if (t_neg > 0 && t_pos > t_neg && i != id) {
	    closest_t = t_neg;
	    sphere_id = i;
	  }
	}
	else {
	  if (t_pos > 0 && t_neg > t_pos && closest_t > t_pos && i != id) {
	    closest_t = t_pos;
	    sphere_id = i;
	  }
	  else if (t_neg > 0 && t_pos > t_neg && closest_t > t_neg && i != id) {
	    closest_t = t_neg;
	    sphere_id = i;
	  }
	}
      }
    }
  }
  if(sphere_id == -1) return 1;
  else return 0;
}











//I = ka*Od + sum(IL * ( kd*Od*(max{0, (N*L)}) + ks*Os*(max{0, (N*H)})^n ))
// computes the color based on the sphere id and light source(s) provided
Pixel * shade_Ray(MTLColorType col[], SphereType sphere[], int sphere_id, int sphereCount, vec3 *eye, RayType *ray, float t, LightType light[], int lightCount, vec3 *viewdir) {
  float shadow = 1.0, il = 1.0;
  //if(sphere_id != 1) printf("%d\n",sphere_id);
  
  //assign the ray/sphere intersection point based on the value t
  vec3* intersection = scaler(ray->dir, t);
  intersection = addVec3(eye, intersection);

  //calcuate the vector N
  vec3* N = subVec3(intersection, sphere[sphere_id].c);
  N = scaler(N, 1 / sphere[sphere_id].r);
  
  float I_r = 0.0, I_g = 0.0, I_b = 0.0;
  I_r += (col->ka * col->odr);
  I_g += (col->ka * col->odg);
  I_b += (col->ka * col->odb);

  //cycle through all light sources assuming at least one light source exists
  //calculate the summation value to be added to the illumination value
  vec3 *L;
  for (int num = 0; num < lightCount; num++) {

    //directional light source
    if (light[num].w == 0.0) {
      L = scaler(light[num].location, -1);
      L = unitVec3(L);
    }
    //point light source
    else{
      L = subVec3(light[num].location, intersection);
      L = unitVec3(L);
    }

    // N * L
    float n_dot_l = dotVec3(N, L);
    if (n_dot_l < 0) n_dot_l = 0.0;

    // Calculating H vector from L and V
    vec3* H_V = scaler(viewdir, -1);
    H_V = addVec3(L, H_V);
    H_V = scaler(H_V, 0.5);
    H_V = unitVec3(H_V);

    // N * H
    float n_dot_h = dotVec3(N, H_V);

    // Taking (N*H) and setting to an exponent
    if (n_dot_h < 0.0) n_dot_h = 0.0;
    else if (n_dot_h != 0.0) {
      int n_abs = abs(col[sphere_id].n);
      
      if (n_abs == 0) n_dot_h = 1.0;
      else if (n_abs > 1) {
	float orig = n_dot_h;
	for (int exp = 2; exp <= n_abs; exp++) {
	  n_dot_h *= orig;
	}
      }
      if (col[sphere_id].n < 0 && n_dot_h != 0.0) n_dot_h = 1 / n_dot_h;
    }

    //Calculate the shadow ray for one specific light source and return the value of S
    shadow = shadow_intersection(light[num].location, light[num].w, intersection, sphere, sphereCount, sphere_id);
    // S * IL
    
    if(lightCount != 1) il = (float) (1.0/lightCount);
    shadow *= il;

    //adding towards the summation as seen in Blinn-Phong illumnation value
    I_r += shadow * ((col->kd * col->odr * n_dot_l) + (col->ks * col->osr * n_dot_h));
    I_g += shadow * ((col->kd * col->odg * n_dot_l) + (col->ks * col->osg * n_dot_h));
    I_b += shadow * ((col->kd * col->odb * n_dot_l) + (col->ks * col->osb * n_dot_h));
    free(L);
  }
  //printf("%f, %f, %f\n", I_r, I_g, I_b);

  //cap values to be between 0 and 1 (inclusive)
  if (I_r > 1) I_r = 1.0;
  else if (I_r < 0) I_r = 0.0;
  if (I_g > 1) I_g = 1.0;
  else if (I_g < 0) I_g = 0.0;
  if (I_b > 1) I_b = 1.0;
  else if (I_b < 0) I_b = 0.0;

  
  //generate color based on ilumination value
  Pixel *color = (Pixel*) calloc(1, sizeof(Pixel));;
  color->r = (int)(255 * I_r);
  color->g = (int)(255 * I_g);
  color->b = (int)(255 * I_b);
  
  return color;
}






// uses ray to check whether there is an intersection between the ray and the array of spheres and returns a color
Pixel * trace_Ray(RayType *ray, vec3 *eye, Pixel * bkgcolor, MTLColorType mtlcolor[], SphereType sphere[], int sphereCount, LightType light[], int lightCount, vec3 *viewdir) {
  float closest_t = 999999999.0;
  int sphere_id = -1;

  //compare each sphere to ray using quadratic formula
  for (int i = 0; i < sphereCount; i++) {
    
    float B = (ray->dir->x * (eye->x - sphere[i].c->x));
    B += (ray->dir->y * (eye->y - sphere[i].c->y));
    B += (ray->dir->z * (eye->z - sphere[i].c->z));
    B *= 2;
    float C = (eye->x - sphere[i].c->x) * (eye->x - sphere[i].c->x);
    C += ((eye->y - sphere[i].c->y) * (eye->y - sphere[i].c->y));
    C += ((eye->z - sphere[i].c->z) * (eye->z - sphere[i].c->z));
    C -= (sphere[i].r * sphere[i].r);
    
    float discriminant = (B * B) - (4 * C);
    if (discriminant > 0) {
      float t_pos = (float) ((-1 * B) + sqrt(discriminant)) / 2;
      float t_neg = (float) ((-1 * B) - sqrt(discriminant)) / 2;

      //assign new closest t value and sphere_id whenevr a intersection occurs t > 0
      if (sphere_id == -1) {
		if (t_pos > 0 && t_pos < t_neg) {
	  closest_t = t_pos;
	  sphere_id = i;
	}
	else if (t_neg > 0 && t_neg < t_pos) {
	  closest_t = t_neg;
	  sphere_id = i;
	}
      }
      else {
	if (t_pos < t_neg && t_pos < closest_t && 0 < t_pos) {
	  closest_t = t_pos;
	  sphere_id = i;
	}
	else if (t_pos > t_neg && t_neg < closest_t && 0 < t_neg) {
	  closest_t = t_neg;
	  sphere_id = i;
	}
      }
    }
  }      
    
  if (sphere_id == -1) {
    return bkgcolor;
  }
  else {
    return shade_Ray(mtlcolor, sphere, sphere_id, sphereCount, eye, ray, closest_t, light, lightCount, viewdir);
  }
}





















//************************************** Main Function **************************************//

int main() {

	using namespace std;
	//************************************** Initialize Main Variables **************************************//
	
	// Information of input scene
	int im_width = -1;
	int im_height = -1;
	vec3 *eye = (vec3*) calloc(1,sizeof(vec3));			// view origin
	vec3 *viewdir = (vec3*) calloc(1,sizeof(vec3));			// Viewing direction
	vec3 *updir = (vec3*) calloc(1,sizeof(vec3));			// Up direction
	float vfov;			// The vertical field of view (degrees)
	Pixel bkgcolor;		// background color




	//************************************** Obtaining and Opening File **************************************//
	
	std::string inputfile;

	// Asks the user whether they wish to use an input file or not
	std::cout << "Please input file name\n";
	std::cin >> inputfile;

	// If the user wishes to input a file with a specific imsize
	std::ifstream myfile;
	myfile.open(inputfile);
	
	if (!myfile.is_open()) {
		std::cout << "Invalid file\n";
		return 1;
	}




	//************************************** Collection of DATA **************************************//

	// Count number of words/strings in file
	int numWords = 0;
	std::string phrase;
	char ch;
	std::string word;
	std::string characters[200];
	while (true) {
	  ch = myfile.get();
	  if (ch == ' ' || ch == '\n'){
	    if(word.length() != 0){
	      characters[numWords] = word;
	      numWords++;
	    }
	    word.clear();
	  }
	  else if(myfile.eof()) {
	    if(word.length() != 0){
	      characters[numWords] = word;
	      numWords++;
	    }
	    word.clear();
	    break;
	  }
	  else if(ch != ' ' && ch != '\n')
	    word += ch;
	}

	int idx = numWords;
	/*
	while (myfile) {
	  std::cout << "Me\n";
	  ch = myfile.get();
	  if(ch != ' ' || ch != '\n'){
	    word += ch;
	  }
	  else{
	    characters[idx] = word;
	    idx++;
	    word.clear();
	  }
	}
	*/
	  
	int sphereCount = 0;
	int materialCount = 0;
	int lightCount = 0;
	//Cycle through the array created from the file's contents and assign values to each input variable
	for (int i = 0; i < idx; i++) {
	  std::string temp = characters[i];
	  if (temp.compare("imsize") == 0) {
	    im_width = std::stoi(characters[i + 1]);
	    im_height = std::stoi(characters[i + 2]);
	    i += 2;
	  }
	  else if (temp.compare("eye") == 0) {
	    eye->x = std::stof(characters[i + 1]);
	    eye->y = std::stof(characters[i + 2]);
	    eye->z = std::stof(characters[i + 3]);
	    i+=3;
	  }
	  else if (temp.compare("viewdir") == 0) {
	    viewdir->x = std::stof(characters[i + 1]);
	    viewdir->y = std::stof(characters[i + 2]);
	    viewdir->z = std::stof(characters[i + 3]);
	    i+=3;
	  }
	  else if (temp.compare("updir") == 0) {
	    updir->x = std::stof(characters[i + 1]);
	    updir->y = std::stof(characters[i + 2]);
	    updir->z = std::stof(characters[i + 3]);
	    i+=3;
	  }
	  else if (temp.compare("vfov") == 0) {
	    vfov = std::stof(characters[i + 1]);
	    i+=1;
	  }
	  else if (temp.compare("bkgcolor") == 0) {
	    bkgcolor.r = (int) (255 * std::stof(characters[i + 1]));
	    bkgcolor.b = (int) (255 * std::stof(characters[i + 2]));
	    bkgcolor.g = (int) (255 * std::stof(characters[i + 3]));
	    i+=3;
	  }
	  else if (temp.compare("sphere") == 0) {
	    sphereCount++;
	    i+=4;
	  }
	  else if (temp.compare("mtlcolor") == 0) {
	    materialCount++;
	    i+=10;
	  }
	  else if (temp.compare("light") == 0){
	    lightCount++;
	    i+=7;
	  }
	  else continue;
	}

	//Cycle through the array and assign input values to arrays (spheres, material colors, and light sources)
	MTLColorType mtlcolor[materialCount+1];	// Material color
	SphereType sphere[sphereCount];	//sphere to be assigned
	LightType light[lightCount];
	int alpha = 0, beta = 0, delta = 0;
	for (int i = 0; i < idx; i++) {
	  std::string temp = characters[i];
	  if (temp.compare("mtlcolor") == 0) {
	    mtlcolor[alpha].odr = std::stof(characters[i + 1]);
	    mtlcolor[alpha].odg = std::stof(characters[i + 2]);
	    mtlcolor[alpha].odb = std::stof(characters[i + 3]);
	    mtlcolor[alpha].osr = std::stof(characters[i + 4]);
	    mtlcolor[alpha].osg = std::stof(characters[i + 5]);
	    mtlcolor[alpha].osb = std::stof(characters[i + 6]);
	    mtlcolor[alpha].ka = std::stof(characters[i + 7]);
	    mtlcolor[alpha].kd = std::stof(characters[i + 8]);
	    mtlcolor[alpha].ks = std::stof(characters[i + 9]);
	    mtlcolor[alpha].n = std::stof(characters[i + 10]);
	    
	    alpha++;
	    i+=10;
	  }
	  else if (temp.compare("sphere") == 0) {
	    sphere[beta].c = (vec3*) calloc(1,sizeof(vec3));
	    sphere[beta].c->x = std::stof(characters[i + 1]);
	    sphere[beta].c->y = std::stof(characters[i + 2]);
	    sphere[beta].c->z = std::stof(characters[i + 3]);
	    sphere[beta].r = std::stof(characters[i + 4]);
	    sphere[beta].m = alpha;
	    beta++;
	    i+=4;
	  }
	  else if(temp.compare("light") == 0){
	    light[delta].location = (vec3*) calloc(1,sizeof(vec3));
	    light[delta].pix = (vec3*) calloc(1,sizeof(vec3));
	    light[delta].location->x = std::stof(characters[i + 1]);
	    light[delta].location->y = std::stof(characters[i + 2]);
	    light[delta].location->z = std::stof(characters[i + 3]);
	    light[delta].w = std::stof(characters[i + 4]);
	    light[delta].pix->x = std::stof(characters[i + 5]);
	    light[delta].pix->y = std::stof(characters[i + 6]);
	    light[delta].pix->z = std::stof(characters[i + 7]);
	    delta++;
	    i+=7;
	  }
	  else continue;
	}
	





	

	//************************************** Variable Error Testing **************************************//
	
	// If width or height are negative then process ends
	if (im_width <= 0 || im_height <= 0) {
	  std::cout << "ERROR: File does not contain proper imsize\n";
	  return 1;
	}

	if (lengthVec3(viewdir) == 0) {
	  std::cout << "ERROR: view direction is invalid\n";
	  return 1;
	}

	if (lengthVec3(updir) == 0) {
	  std::cout << "ERROR: up direction is invalid\n";
	  return 1;
	}


	



	
	

	//************************************** Viewing Window Calculation **************************************//
	
	// Need to use the updir, viewdir, and eye to calculate the four corners of the viewing window
	
	// Defining the Camera coordinate system
	vec3 *v, *u, *w;
	w = scaler(viewdir, -1);

	vec3 *viewTemp = crossProdVec3(viewdir,updir);
	if (lengthVec3(viewTemp) == 0) {
		std::cout << "ERROR: viewdir == updir. Please redefine\n";
		return 1;
	}
	float viewTempU = lengthVec3(viewTemp);
	u = scaler(viewTemp, (1/viewTempU));

	vec3 *viewTemp1 = crossProdVec3(u,viewdir);
	viewTempU = lengthVec3(viewTemp1);
	v = scaler(viewTemp1, (1/viewTempU));
	free(viewTemp1);
	
	free(viewTemp);
	


	
	//Defining the four corners of the viewing window
	vec3 *ul, *ur, *ll, *lr; 

	float d = 10; // arbitrary value used to calculate the viewing window

	float height = tan((vfov/2)) * (2 * d);
	float aspect_ratio = (float)im_width / (float)im_height;
	float width = height * aspect_ratio;

	viewTempU = lengthVec3(viewdir);
	vec3 *cornerTemp = scaler(viewdir, (1/viewTempU));
	//value1 = eye     value2 = cornerTemp1    value3 = cornerTemp2     value4 = cornerTemp3
	vec3 *cornerTemp1 = scaler(cornerTemp,d);
	vec3 *cornerTemp2 = scaler(u, width/2);
	vec3 *cornerTemp3 = scaler(v, height/2);

	//eye + scaler(cornerTemp, d) + scaler(u, width/2) + scaler(v, height/2)
	
	//ul = addVec3(subVec3(addVec3(eye, scaler(unitVec3(viewdir), d)), scaler(u, width / 2)), scaler(v, height / 2));
	ul = addVec3(eye,cornerTemp1);
	ul = subVec3(ul,cornerTemp2);
	ul = addVec3(ul,cornerTemp3);
	
	//ur = addVec3(addVec3(addVec3(eye, scaler(unitVec3(viewdir), d)), scaler(u, width / 2)), scaler(v, height / 2));
	ur = addVec3(eye,cornerTemp1);
	ur = addVec3(ur,cornerTemp2);
	ur = addVec3(ur,cornerTemp3);
	
	//ll = subVec3(subVec3(addVec3(eye, scaler(unitVec3(viewdir), d)), scaler(u, width / 2)), scaler(v, height / 2));
	ll = addVec3(eye,cornerTemp1);
	ll = subVec3(ll,cornerTemp2);
	ll = subVec3(ll,cornerTemp3);
	
	//lr = subVec3(addVec3(addVec3(eye, scaler(unitVec3(viewdir), d)), scaler(u, width / 2)), scaler(v, height / 2));
	lr = addVec3(eye,cornerTemp1);
	lr = addVec3(lr,cornerTemp2);
	lr = subVec3(lr,cornerTemp3);

	free(cornerTemp);
	free(cornerTemp1);
	free(cornerTemp2);
	free(cornerTemp3);
	


	


	
	//************************************** Mapping between Pixels and Viewing Window **************************************//
	
	vec3 *deltaH, *deltaV, *deltaC_H, *deltaC_V;
	vec3 *tempDelta = subVec3(ur,ul);
	vec3 *tempDelta1 = subVec3(ll,ul);

	float iwid = (float) 1/im_width;
	float ihei = (float) 1/im_height;

	deltaH = scaler(tempDelta, iwid);
	deltaV = scaler(tempDelta1, ihei);

	deltaC_H = scaler(tempDelta, (0.5*iwid));
	deltaC_V = scaler(tempDelta1,(0.5*ihei));

	free(tempDelta);
	free(tempDelta1);
	
    


	

	//************************************** Representation of each Ray and Color Assignment **************************************//

	/*
		ul + (i)*deltaH + (j)*deltaV + deltaC_H + deltaC_V
		i = 0 ... (im_width - 1)
		j = 0 ... (im_height - 1)
	*/

	RayType *ray = (RayType*) calloc(1, sizeof(RayType));
	Pixel *img = (Pixel*) calloc(im_width*im_height, sizeof(Pixel));

	for (uint32_t i = 0; i < im_width; i++) {
	  for (uint32_t j = 0; j < im_height; j++) {
	    vec3 *temp1 = scaler(deltaH, i);
	    vec3 *temp2 = scaler(deltaV, j);

	    //RayType *ray = (RayType*) calloc(1, sizeof(RayType));
	    ray->p = ul;
	    ray->p = addVec3(ray->p, temp1);
	    ray->p = addVec3(ray->p, temp2);
	    ray->p = addVec3(ray->p, deltaC_H);
	    ray->p = addVec3(ray->p, deltaC_V);
	    
	    vec3 *temp4 = subVec3(ray->p,eye);
	    float tempValue = lengthVec3(temp4);
	    tempValue = (float) (1/tempValue);
	    ray->dir = scaler(temp4, tempValue);

		
	    
	    //printf("%d, %d\n", i, j);
	    
	    Pixel *temp = trace_Ray(ray, eye, &bkgcolor, mtlcolor, sphere, sphereCount, light, lightCount, viewdir);
	    img[i*im_width+j].r = temp->r;
	    img[i*im_width+j].g = temp->g;
	    img[i*im_width+j].b = temp->b;

	    if(temp != &bkgcolor){
	      free(temp);
	    }
	    
	    free(temp1);
	    free(temp2);
	  }
	}
	free(ray);
        
	

	//************************************** Generation of PPM File **************************************//
	
	// Finds the suffix of the input file and generates the output filename with .ppm
	int val = inputfile.find(".");
	//int len = inputfile.length() - val;
	std::string outputfile = inputfile.substr(0, val) + ".ppm";

	std::ofstream output_stream(outputfile, std::ios::out);

	output_stream << "P3\n"
		      << im_width << "\n"
		      << im_height << "\n"
		      << 255 << "\n";
	for (uint32_t x = 0; x < im_width; x++){
	  for (uint32_t y = 0; y < im_height; y++){
	    output_stream << (unsigned int)img[x*im_width + y].r << " "
			  << (unsigned int)img[x*im_width + y].g << " "
			  << (unsigned int)img[x*im_width + y].b << "\n";
   
	    if((x+1) == im_width && (y+1)==im_height) break;
	  }
	  
	}

	output_stream.close();
	exit(0);
	return 0;
}

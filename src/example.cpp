#include "camera.hpp"

#include <iostream>
#include <bits/stdc++.h>
#include "glm/ext.hpp"
#include "glm/gtx/string_cast.hpp"

using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;
using namespace std;

GL::Rasterizer r;
GL::ShaderProgram program;

const int n = 20;
const int nv = n*n;
const int nt = 2*(n-1)*(n-1);
vec3 vertices[nv];
vec3 normals[nv];
ivec3 triangles[nt];

const int n1 = 20;
const int nv1 = n1*n1;
const int nt1 = 2*(n1-1)*(n1-1);
vec3 vertices1[nv1];
vec3 normals1[nv1];
ivec3 triangles1[nt1];


GL::Object object;
GL::AttribBuf vertexBuf, normalBuf;

// GL::Rasterizer r1;
// GL::ShaderProgram program1;

GL::Object object1;
GL::AttribBuf vertexBuf1, normalBuf1;

//vec3 vertices2[];
//vec3 normals2[];
//ivec3 triangles1[nt1];
// #define PI glm::pi <GLfloat> ()
// #define stacks 30
// #define slices 30
// #define radius 2
// int numberOfPoints = 0;
// int numberOfIndexes = 0;
// vec3 vertices2[(stacks + 1) * slices * 3];
// ivec3 indices2[stacks * slices * 10];
// int dx = 0, dy = 0, dz = 0;

GL::Object object2;
GL::AttribBuf vertexBuf2, normalBuf2;

CameraControl camCtl;

vector<vec3> position(n*n);
vector<vec3> permanent(n*n);
vector<vec3> velocity(n*n);

vector<vec3> position1(n1*n1);
vector<vec3> permanent1(n1*n1);
vector<vec3> velocity1(n1*n1);

int k_struct = 500;
int k_shear = 100;
int k_flex = 10;


vector<int> is_fixed(n*n,0);

float k_struct_dist = 1/n;
float k_shear_dist = 1/n*sqrt(2);
float k_flex_dist = 2/n;


float damp = 0.3;

float mass = 0.1;

void initializeScene() {
	object = r.createObject();

	// cout << "nv = " << 	nv << endl;
	for(int i = 0 ; i<nv ; i++){
		int toy = i/n;
		float tt = (float) toy;
		vertices[i] = position[i] = vec3(((float) toy)/((float) n) , 0.0 ,(float)(i%n)/n);
		// cout << i/(n*n) << " " << (i%n)/n << "\n";
	}

	for(int i = 0 ; i<nv ; i++){
		permanent[i] = position[i];
	}


	vertexBuf = r.createVertexAttribs(object, 0, nv, vertices);

	for(int i = 0; i<nv ; i++){
		normals[i] = vec3(0,0,1);
	}

	for(int i = 0; i<nv ; i++){
		velocity[i] = vec3(0.0, 0.0, 0.0);
	}

	// for(int i = 0; i < n ; i++){
	// 	is_fixed[i*n] = 1;
	// }
	is_fixed[0] = 1;
	is_fixed[(n-1)*n] = 1;

	normalBuf = r.createVertexAttribs(object, 1, nv, normals);

	int cnt = 0;
	for(int i = 0; i<n-1 ; i++){
		for(int j = 0; j<n -1 ; j++){
			triangles[cnt] = ivec3(i*n + j, i*n+j +1, (i+1)*n+j);
			cnt++;
			triangles[cnt] = ivec3(i*n + j +1, (i+1)*n + j + 1, (i+1)*n+j);
			cnt++;
		}
	}
	r.createTriangleIndices(object, nt, triangles); 

	object1 = r.createObject();

	//cout << "nv = " << 	nv1 << endl;
	for(int i = 0 ; i<nv1 ; i++){
		int toy = i/n1;
		float tt = (float) toy;
		vertices1[i] = position1[i] = vec3(((float) toy)/((float) n1), -0.9,(float)(i%n1)/n1-0.25);
		//cout << i/(n1*n1) << " " << (i%n1)/n1 << "\n";
	}

	for(int i = 0 ; i<nv1 ; i++){
		permanent1[i] = position1[i];
	}


	vertexBuf1 = r.createVertexAttribs(object1, 0, nv1, vertices1);

	for(int i = 0; i<nv1 ; i++){
		normals1[i] = vec3(0,0,1);
	}

	for(int i = 0; i<nv1 ; i++){
		velocity1[i] = vec3(0.0, 0.0, 0.0);
	}

	normalBuf1 = r.createVertexAttribs(object1, 1, nv1, normals1);

	int cnt1 = 0;
	for(int i = 0; i<n1-1 ; i++){
		for(int j = 0; j<n1 -1 ; j++){
			triangles1[cnt1] = ivec3(i*n1 + j, i*n1+j +1, (i+1)*n1+j);
			cnt1++;
			triangles1[cnt1] = ivec3(i*n1 + j +1, (i+1)*n1 + j + 1, (i+1)*n1+j);
			cnt1++;
		}
	}
	r.createTriangleIndices(object1, nt1, triangles1); 

	// object2 = r.createObject();

	// for (int i = 0; i <= slices; i++)
    // {
    //     float phi = i * (glm::pi<float>() / slices) ;

    //     for (int j = 0; j < stacks; j++)
    //     {
    //         float theta = j * (glm::pi <float>() * 2 / stacks);

    //         float x = cosf(theta) * sinf(phi);
    //         float y = cosf(phi);
    //         float z = sinf(theta) * sinf(phi);

    //         vertices[numberOfPoints++] = x * radius;
    //         vertices[numberOfPoints++] = y * radius;
    //         vertices[numberOfPoints++] = z * radius;
    //     }
    // }
    // for (int i = 0; i < numberOfPoints / 3 - stacks; i++)
    // {
    //     if ((i + 1) % stacks == 0)
    //     {
    //         indices[numberOfIndexes++] = i;
    //         indices[numberOfIndexes++] = i - stacks + 1;
    //         indices[numberOfIndexes++] = i + stacks;

    //         indices[numberOfIndexes++] = i - stacks + 1;
    //         indices[numberOfIndexes++] = i + stacks;
    //         if (i + 1 == numberOfPoints / 3)
    //             indices[numberOfIndexes++] = numberOfPoints - stacks;
    //         else
    //             indices[numberOfIndexes++] = i + 1;
    //     }
    //     else
    //     {
    //         indices[numberOfIndexes++] = i;
    //         indices[numberOfIndexes++] = i + 1;
    //         indices[numberOfIndexes++] = i + stacks;

    //         indices[numberOfIndexes++] = i + 1;
    //         indices[numberOfIndexes++] = i + stacks;
    //         indices[numberOfIndexes++] = i + stacks + 1;
    //     }
    // }
}


//collision check with plane
// bool collisionchecK(){
// 	for(int i =0;i<nv;i++){
// 		if(position[i].y<=-1.0){
// 			//collision true
// 			return true;
// 		}
// 	}

// }
//perfectly placstic collision with plane
void collisionchecK(){
	for(int i =0;i<nv;i++){
		if(position[i].y<=-0.9){
			position[i].y=-0.9;
			//velocity[i]=vec3(0.0,0.0,0.0);
			velocity[i].y=0.0;
			velocity[i].z=(velocity[i].z)/(1.11);
			//v=u.(sin-(e+1)cos)
			
		}
	}

}
void updateScene(float t) {

	for(int i = 0; i<nv ; i++){
		if(is_fixed[i] == 1 ){
			continue;
		}
		vec3 net_Force = vec3(0.0, -mass*10, 0.0);
		// vec3 net_Force = vec3(0.0, 0.0, 0.0);

		int posx = i/n;
		int posy = i%n;

		vec3 net_Force2 = vec3(0.0,0.0,0.0);
		vec3 right = vec3(0.0,0.0,0.0); vec3 left = vec3(0.0,0.0,0.0); 
		vec3 up = vec3(0.0,0.0,0.0); vec3 down = vec3(0.0,0.0,0.0);

		// STRUCTURAL SPRINGS

		if(posy<n-1){
			right = position[(posx)*n + posy + 1];
			right = right - position[i];
			float l = sqrt(right.x*right.x + right.y*right.y + right.z*right.z);
			vec3 temp = permanent[(posx)*n + posy + 1] - permanent[i];
			float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
			l = l - ll;
			if(l < 1e-5){
				l = 0;
			}
			l = l*k_struct;
			right = normalize(right);
			net_Force2 = net_Force2 + vec3(right.x*l, right.y*l, right.z*l);
		}
		if(posy > 0){
			left = position[(posx)*n + posy -1];
			left = left - position[i];
			float l = sqrt(left.x*left.x + left.y*left.y + left.z*left.z);
			vec3 temp = permanent[(posx)*n + posy - 1] - permanent[i];
			float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
			l = l - ll;
			if(l < 1e-5){
				l = 0;
			}
			l = l*k_struct;
			left = normalize(left);
			net_Force2 = net_Force2 + vec3(left.x*l, left.y*l, left.z*l);
		}
		if(posx>0){
			up = position[(posx-1)*n + posy];
			up = up - position[i];
			float l = sqrt(up.x*up.x + up.y*up.y + up.z*up.z);
			vec3 temp = permanent[(posx-1)*n + posy] - permanent[i];
			float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
			l = l - ll;
			if(l < 1e-5){
				l = 0;
			}
			l = l*k_struct;
			up = normalize(up);
			net_Force2 = net_Force2 + vec3(up.x*l, up.y*l, up.z*l);
		}
		if(posx<n-1){
			down = position[(posx+1)*n + posy];
			down = down - position[i];
			float l = sqrt(down.x*down.x + down.y*down.y + down.z*down.z);
			vec3 temp = permanent[(posx+1)*n + posy] - permanent[i];
			float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
			l = l - ll;
			if(l < 1e-5){
				l = 0;
			}
			l = l*k_struct;
			// l = 0;
			down = normalize(down);
			net_Force2 = net_Force2 + vec3(down.x*l, down.y*l, down.z*l);
		}

		

		// SHEAR SPRINGS

		if(posy<n-1 && posx < n-1){
			right = position[(posx+1)*n + posy + 1];
			right = right - position[i];
			float l = sqrt(right.x*right.x + right.y*right.y + right.z*right.z);
			vec3 temp = permanent[(posx+1)*n + posy+1] - permanent[i];
			float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
			l = l - ll;
			if(l < 1e-5){
				l = 0;
			}
			l = l*k_shear;
			right = normalize(right);
			net_Force2 = net_Force2 + vec3(right.x*l, right.y*l, right.z*l);
		}
		if(posy > 0 && posx < n-1){
			left = position[(posx+1)*n + posy -1];
			left = left - position[i];
			float l = sqrt(left.x*left.x + left.y*left.y + left.z*left.z);
			vec3 temp = permanent[(posx+1)*n + posy+1] - permanent[i];
			float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
			l = l - ll;
			if(l < 1e-5){
				l = 0;
			}
			l = l*k_shear;
			// l = 0;
			left = normalize(left);
			net_Force2 = net_Force2 + vec3(left.x*l, left.y*l, left.z*l);
		}
		if(posx>0 && posy < n-1){
			up = position[(posx-1)*n + posy + 1];
			up = up - position[i];
			float l = sqrt(up.x*up.x + up.y*up.y + up.z*up.z);
			vec3 temp = permanent[(posx-1)*n + posy+1] - permanent[i];						
			float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
			l = l - ll;
			if(l < 1e-5){
				l = 0;
			}
			l = l*k_shear;
			up = normalize(up);
			net_Force2 = net_Force2 + vec3(up.x*l, up.y*l, up.z*l);
		}
		if(posx>0 && posy >0){
			down = position[(posx-1)*n + posy-1];
			down = down - position[i];
			float l = sqrt(down.x*down.x + down.y*down.y + down.z*down.z);
			vec3 temp = permanent[(posx-1)*n + posy-1] - permanent[i];
			float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
			l = l - ll;
			if(l < 1e-5){
				l = 0;
			}
			l = l*k_shear;
			down = normalize(down);
			net_Force2 = net_Force2 + vec3(down.x*l, down.y*l, down.z*l);
		}

		// BENDING SPRINGS

		if(posy<n-2){
			right = position[(posx)*n + posy + 2];
			right = right - position[i];
			float l = sqrt(right.x*right.x + right.y*right.y + right.z*right.z);
			vec3 temp = permanent[(posx)*n + posy+2] - permanent[i];
			float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
			l = l - ll;
			if(l < 1e-5){
				l = 0;
			}
			l = l*k_flex;
			right = normalize(right);
			net_Force2 = net_Force2 + vec3(right.x*l, right.y*l, right.z*l);
		}
		if(posy > 1){
			left = position[(posx)*n + posy -2];
			left = left - position[i];
			float l = sqrt(left.x*left.x + left.y*left.y + left.z*left.z);
			vec3 temp = permanent[(posx)*n + posy-2] - permanent[i];
			float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
			l = l - ll;
			if(l < 1e-5){
				l = 0;
			}
			l = l*k_flex;
			left = normalize(left);
			net_Force2 = net_Force2 + vec3(left.x*l, left.y*l, left.z*l);
		}
		if(posx>1){
			up = position[(posx-2)*n + posy];
			up = up - position[i];
			float l = sqrt(up.x*up.x + up.y*up.y + up.z*up.z);
			vec3 temp = permanent[(posx-2)*n + posy] - permanent[i];
			float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
			l = l - ll;
			if(l < 1e-5){
				l = 0;
			}
			l = l*k_flex;
			up = normalize(up);
			net_Force2 = net_Force2 + vec3(up.x*l, up.y*l, up.z*l);
		}
		if(posx<n-2){
			down = position[(posx+2)*n + posy];
			down = down - position[i];
			float l = sqrt(down.x*down.x + down.y*down.y + down.z*down.z);
			vec3 temp = permanent[(posx+2)*n + posy] - permanent[i];
			float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
			l = l - ll;
			if(l < 1e-5){
				l = 0;
			}
			l = l*k_flex;
			down = normalize(down);
			net_Force2 = net_Force2 + vec3(down.x*l, down.y*l, down.z*l);
		}

		net_Force += net_Force2;
		net_Force += vec3(velocity[i].x * (-1*damp), velocity[i].y *(-1*damp), velocity[i].z *(-1*damp));
		// std::cout << glm::to_string(net_Force2) << std::endl;
		vec3 accn = vec3(net_Force.x/mass, net_Force.y /mass, net_Force.z /mass); 
		velocity[i] += vec3(accn.x *t, accn.y *t, accn.z *t);
		position[i] += vec3(velocity[i].x *t, velocity[i].y *t, velocity[i].z *t);
	}
	//CONSTRAINT PROJECTION

    // int loop_iterations = 10;
    // while(loop_iterations--){
    //     for(int i = 0; i<nv ; i++){
    //         int posx = i/n;
    //         int posy = i%n;
    //         vec3 right, left, up , down;

        
    //         // I'm safe from right side
    //         if(posy<n-1){
    //             right = position[(posx)*n + posy + 1];
    //             right = right - position[i];
    //             float l = sqrt(right.x*right.x + right.y*right.y + right.z*right.z);
    //             vec3 temp = permanent[(posx)*n + posy + 1] - permanent[i];
    //             float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
    //             l = l - ll;
    //             if(l < 1e-5){
    //                 l = 0;
    //             }
    //             // l=0;
    //             right = normalize(right);
    //             vec3 p1_change = vec3(-0.5*l*right.x,-0.5*l*right.y,-0.5*l*right.z);

    //             if(is_fixed[i]==1){
    //                 position[(posx)*n + posy + 1] = position[(posx)*n + posy + 1] + vec3(p1_change.x*2, p1_change.y*2, p1_change.z*2 );
    //             }
    //             else if(is_fixed[(posx)*n + posy+1] == 1){
    //                 position[i] = position[i] - vec3(p1_change.x*2, p1_change.y*2, p1_change.z*2 );
    //             }
    //             else{
    //                 position[i] = position[i] - p1_change;
    //                 position[(posx)*n + posy + 1] = position[(posx)*n + posy + 1] + p1_change;
    //             }
    //         }
    //         // // I'm safe from left side
    //         // if(posy > 0){
    //         //  left = position[(posx)*n + posy -1];
    //         //  left = left - position[i];
    //         //  float l = sqrt(left.x*left.x + left.y*left.y + left.z*left.z);
    //         //  vec3 temp = permanent[(posx)*n + posy - 1] - permanent[i];
    //         //  float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
    //         //  l = l - ll;
    //         //  if(l < 1e-5){
    //         //      l = 0;
    //         //  }
    //         //  // l = 0;
    //         //  // l = l*k_struct;
    //         //  left = normalize(left);
    //         //  vec3 p1_change = vec3(-0.5*l*left.x,-0.5*l*left.y,-0.5*l*left.z);

                
    //         //  if(is_fixed[i]==1){
    //         //      position[(posx)*n + posy - 1] = position[(posx)*n + posy - 1] + vec3(p1_change.x*2, p1_change.y*2, p1_change.z*2 );
    //         //  }
    //         //  else{
    //         //      position[i] = position[i] - p1_change;
    //         //      position[(posx)*n + posy - 1] = position[(posx)*n + posy - 1] + p1_change;
    //         //  }
    //         //  // net_Force2 = net_Force2 + vec3(left.x*l, left.y*l, left.z*l);
    //         // }
    //         // Safe from above
    //         if(posx>0){
    //             up = position[(posx-1)*n + posy];
    //             up = up - position[i];
    //             float l = sqrt(up.x*up.x + up.y*up.y + up.z*up.z);
    //             vec3 temp = permanent[(posx-1)*n + posy] - permanent[i];
    //             float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
    //             l = l - ll;
    //             if(l < 1e-5){
    //                 l = 0;
    //             }
    //             // l = l*k_struct;
    //             up = normalize(up);
    //             vec3 p1_change = vec3(-0.5*l*up.x,-0.5*l*up.y,-0.5*l*up.z);

    //             if(is_fixed[i]==1){
    //                 position[(posx-1)*n + posy] = position[(posx-1)*n + posy] + vec3(p1_change.x*2, p1_change.y*2, p1_change.z*2 );
    //             }
    //             else if(is_fixed[(posx-1)*n + posy] == 1){
    //                 position[i] = position[i] - vec3(p1_change.x*2, p1_change.y*2, p1_change.z*2 );
    //             }
    //             else{
    //                 position[i] = position[i] - p1_change;
    //                 position[(posx-1)*n + posy] = position[(posx-1)*n + posy] + p1_change;
    //             }
    //         }
    //         // // Safe from below
    //         // if(posx<n-1){
    //         //  down = position[(posx+1)*n + posy];
    //         //  down = down - position[i];
    //         //  float l = sqrt(down.x*down.x + down.y*down.y + down.z*down.z);
    //         //  vec3 temp = permanent[(posx+1)*n + posy] - permanent[i];
    //         //  float ll = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
    //         //  l = l - ll;
    //         //  if(l < 1e-5){
    //         //      l = 0;
    //         //  }
    //         //  // l = l*k_struct;
    //         //  // l = 0;
    //         //  down = normalize(down);
    //         //  vec3 p1_change = vec3(-0.5*l*down.x,-0.5*l*down.y,-0.5*l*down.z);
                
    //         //  if(is_fixed[i]==1){
    //         //      position[(posx+1)*n + posy] = position[(posx+1)*n + posy] + vec3(p1_change.x*2, p1_change.y*2, p1_change.z*2 );
    //         //  }
    //         //  else{
    //         //      position[i] = position[i] - p1_change;
    //         //      position[(posx+1)*n + posy] = position[(posx+1)*n + posy] + p1_change;
    //         //  }
    //         //  // net_Force2 = net_Force2 + vec3(down.x*l, down.y*l, down.z*l);
    //         // }

    //     }
    // }


	collisionchecK();

	for(int i = 0 ; i<nv ; i++){
		vertices[i] = position[i];
	}
	

	vertexBuf = r.createVertexAttribs(object, 0, nv, vertices);

	for(int i = 0; i<nv ; i++){
        int posx = i/n;
        int posy = i%n;
        if(posx<n-1 && posy<n-1){
            normals[i] = normalize(cross(vertices[(posx+1)*n+posy] - vertices[i], vertices[i] - vertices[posx*n+posy+1]));
        }
        else if(posx<n-1 && posy==n-1){
            normals[i] = normalize(cross(vertices[(posx+1)*n+posy] - vertices[i], vertices[posx*n+posy-1] - vertices[i]));
        }
        else if(posx==n-1 && posy<n-1){
            normals[i] = normalize(cross( vertices[(posx-1)*n+posy] - vertices[i], vertices[posx*n+posy+1] - vertices[i]));
        }
        else{
            normals[i] = normalize(cross( vertices[(posx-1)*n+posy] - vertices[i],  vertices[i] - vertices[posx*n+posy-1] ));
        }
    }
    r.updateVertexAttribs(normalBuf, nv, normals);
}


int main() {
	int width = 640, height = 480;
	if (!r.initialize("Animation", width, height)) {
		return EXIT_FAILURE;
	}
	camCtl.initialize(width, height);
	camCtl.camera.setCameraView(vec3(0.5, -0.5, 1.5), vec3(0.5, -0.5, 0.0), vec3(0.0, 1.0, 0.0));
	program = r.createShaderProgram(
		r.vsPhongShading(),
		r.fsPhongShading()
	);
	// program1 = r1.createShaderProgram(
	// 	r1.vsPhongShading(),
	// 	r1.fsPhongShading()
	// );

	initializeScene();

	float prev = 0.0;
	while (!r.shouldQuit()) {
        //float t = SDL_GetTicks64()*1e-3;
		updateScene(0.001);
		// cout << t << ";

		//prev = t;

		camCtl.update();
		Camera &camera = camCtl.camera;

		r.clear(vec4(1.0, 1.0, 1.0, 1.0));
		r.enableDepthTest();
		r.useShaderProgram(program);

		r.setUniform(program, "model", glm::mat4(1.0));
		r.setUniform(program, "view", camera.getViewMatrix());
		r.setUniform(program, "projection", camera.getProjectionMatrix());
		r.setUniform(program, "lightPos", camera.position);
		r.setUniform(program, "viewPos", camera.position);
		r.setUniform(program, "lightColor", vec3(1.0f, 1.0f, 1.0f));

		r.setupFilledFaces();
		r.setUniform(program, "objectColor", vec3(1.0f, 0.5f, 0.0f));
		r.drawObject(object1);
		r.drawObject(object);
		r.setupWireFrame();
		r.setUniform(program, "objectColor", vec3(0.0f, 0.0f, 0.0f));
		r.drawObject(object1);
		r.drawObject(object);
		r.show();
		// r.setupFilledFaces();
		// r.setUniform(program, "objectColor", vec3(1.0f, 1.0f, 0.0f));
		// r.drawObject(object1);

		// r.setupWireFrame();
		// r.setUniform(program, "objectColor", vec3(0.0f, 0.0f, 0.0f));
		// r.drawObject(object1);
		// r1.clear(vec4(1.0, 1.0, 1.0, 1.0));
		// r1.enableDepthTest();
		// r1.useShaderProgram(program1);

		// r1.setUniform(program1, "model", glm::mat4(1.0));
		// r1.setUniform(program1, "view", camera.getViewMatrix());
		// r1.setUniform(program1, "projection", camera.getProjectionMatrix());
		// r1.setUniform(program1, "lightPos", camera.position);
		// r1.setUniform(program1, "viewPos", camera.position);
		// r1.setUniform(program1, "lightColor", vec3(1.0f, 1.0f, 1.0f));

		// r1.setupFilledFaces();
		// r1.setUniform(program1, "objectColor", vec3(1.0f, 0.5f, 0.0f));
		// r1.drawObject(object1);
		// r1.setupWireFrame();
		// r1.setUniform(program1, "objectColor", vec3(0.0f, 0.0f, 0.0f));
		// r1.drawObject(object1);

		// r1.show();
	}
}
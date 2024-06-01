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

const int n = 10;
const int nv = n*n;
const int nt = 2*(n-1)*(n-1);
vec3 vertices[nv];
vec3 normals[nv];

// float lats = 40;
// float longs = 40;
// float m_vao = 0;
// float m_vboVertex = 0;
// float m_vboIndex = 0;
// vec3 vertices1[lats*longs];
// vec3 normals1[lats*longs];

ivec3 triangles[nt];

GL::Object object;
GL::AttribBuf vertexBuf, normalBuf;

// GL::Object object1;
// GL::AttribBuf vertexBuf1, normalBuf1;

CameraControl camCtl;

vector<vec3> position(n*n);
vector<vec3> permanent(n*n);
vector<vec3> velocity(n*n);

int k_struct = 50;
int k_shear = 0;
int k_flex = 0;


vector<int> is_fixed(n*n,0);

float k_struct_dist = 1/n;
float k_shear_dist = 1/n*sqrt(2);
float k_flex_dist = 2/n;


float damp = 0.4;

float mass = 0.1;

void initializeScene() {
	object = r.createObject();

	// cout << "nv = " << 	nv << endl;
	for(int i = 0 ; i<nv ; i++){
		int toy = i/n;
		float tt = (float) toy;
		vertices[i] = position[i] = vec3(((float) toy)/((float) n) ,0.0 ,(float)(i%n)/n);
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

	// object1=r.createObject();

	// float radius = 2.0;
	// float PI = 3.14;
	// float sectorCount = 36;
	// float stackCount = 18;
	// float x, y, z, xy;                              // vertex position
	// float nx, ny, nz, lengthInv = 1.0f / radius;    // vertex normal
	// float s, t;                                     // vertex texCoord

	// float sectorStep = 2 * PI / sectorCount;
	// float stackStep = PI / stackCount;
	// float sectorAngle, stackAngle;
	// for(int i = 0; i <= stackCount; ++i){
	// 	stackAngle = PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
	// 	xy = radius * cosf(stackAngle);             // r * cos(u)
	// 	z = radius * sinf(stackAngle);              // r * sin(u)

	// 	// add (sectorCount+1) vertices per stack
	// 	// first and last vertices have same position and normal, but different tex coords
	// 	for(int j = 0; j <= sectorCount; ++j){
	// 		sectorAngle = j * sectorStep;           // starting from 0 to 2pi

	// 		// vertex position (x, y, z)
	// 		x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
	// 		y = xy * sinf(sectorAngle);             // r * cos(u) * sin(v)
	// 		vertices1.push_back(x);
	// 		vertices1.push_back(y);
	// 		vertices1.push_back(z);

	// 		vertexBuf1 = r.createVertexAttribs(object1,3,vertices1.size(), vertices1);

	// 		// normalized vertex normal (nx, ny, nz)
	// 		nx = x * lengthInv;
	// 		ny = y * lengthInv;
	// 		nz = z * lengthInv;
	// 		normals1.push_back(nx);
	// 		normals1.push_back(ny);
	// 		normals1.push_back(nz);

	// 		normalBuf1 = r.createVertexAttribs(object1,3, radius, normals1);
	// 		// vertex tex coord (s, t) range between [0, 1]
	// 		// s = (float)j / sectorCount;
	// 		// t = (float)i / stackCount;
	// 		// texCoords.push_back(s);
	// 		// texCoords.push_back(t);
	// 	}
	// }

// 	object1=r.createObject();

// 	float pos[]={
        
//         3.0f,0.0f,2.0f,
//         3.0f,0.0f,4.0f,
//         5.0f,0.0f,4.0f,
//         5.0f,0.0f,2.0f,
        
        
//         3.0f,2.0f,2.0f,
//         3.0f,2.0f,4.0f,
//         5.0f,2.0f,4.0f,
//         5.0f,2.0f,2.0f
        
//     };
    
    
// //     float color[]={
        
// //         1.0f,0.0f,1.0f,
// //         0.6f,0.2f,0.3f,
// //         0.9f,0.3f,0.0f,
// //         0.3f,0.9f,0.2f,
        
        
// //         1.0f,0.0f,1.0f,
// //         0.6f,0.2f,0.3f,
// //         0.9f,0.3f,0.0f,
// //         0.3f,0.9f,0.2f
        
// //     };
    
// // 	//    float normal[]={
// // 	//
// // 	//
// // 	//    };
    
// //    unsigned int indices[]{
// //         1,0,2,
// //         2,0,3,
// //         1,2,5,
// //         2,6,5,
// //         3,7,6,
// //         3,6,2,
// //         6,4,5,
// //         6,7,4,
// //         1,5,4,
// //         1,4,0,
// //         0,4,7,
// //         3,0,7
// //     };

// 	vertexBuf1 = r.createVertexAttribs(object, 0,8, pos);

	// glBufferData(GL_ARRAY_BUFFER,sizeof(pos),&pos[0],GL_STATIC_DRAW);
    // glEnableVertexAttribArray( 0 );
    // glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

}

void updateScene(float t) {

	for(int i = 0; i<nv ; i++){
		if(is_fixed[i] == 1 ){
			continue;
		}
		//vec3 net_Force = vec3(0.0, -mass*10, 0.0);
		vec3 net_Force = vec3(0.0, 0.0, 0.0);

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
			l = l*k_struct ;
			// l = 0;
			down = normalize(down);
			net_Force2 = net_Force2 + vec3(down.x*l, down.y*l, down.z*l);
		}

		// ----------CONSTRAINT PROJECTION --------------

		

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
	for(int i = 0 ; i<nv ; i++){
		vertices[i] = position[i];
	}
	vertexBuf = r.createVertexAttribs(object, 0, nv, vertices);
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

	initializeScene();

	float prev = 0.0;
	while (!r.shouldQuit()) {
        // float t = SDL_GetTicks64()*1e-3;

		updateScene(0.001);
		// cout << t << ";

		// prev = t;

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
		r.drawObject(object);

		r.setupWireFrame();
		r.setUniform(program, "objectColor", vec3(0.0f, 0.0f, 0.0f));
		r.drawObject(object);

		r.show();
	}
}
#include <Windows.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/freeglut.h>

#define GLFW_INCLUDE_GLU
#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <vector>

#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>

using namespace glm;

// -------------------------------------------------
// Global Variables
// -------------------------------------------------
int Width = 512;
int Height = 512;
std::vector<float> OutputImage;
// -------------------------------------------------

class Ray {
public:
	vec3 origin;
	vec3 direction;

	Ray(const vec3& o, const vec3& d) : origin(o), direction(normalize(d)) {}
};

class Camera {
public:
	vec3 position;
	vec3 u, v, w;
	float l, r, b, t, d;

	Camera(const vec3& pos, const vec3& lookAt, const vec3& up,
		float left, float right, float bottom, float top, float distance)
		: position(pos), l(left), r(right), b(bottom), t(top), d(distance) {

		w = normalize(position - lookAt);
		u = normalize(cross(up, w));
		v = cross(w, u);
	}

	Ray generateRay(float i, float j) const {
		float u_screen = l + (r - l) * (i / Width);
		float v_screen = b + (t - b) * (j / Height);
		vec3 pixelPos = position + (u_screen * u) + (v_screen * v) - (d * w);
		return Ray(position, pixelPos - position);
	}
};

class Surface {
public:
	virtual bool intersect(const Ray& ray, const float& tMin, const float& tMax, float& d) const = 0;
};

class Sphere : public Surface {
public:
	vec3 center;
	float radius;

	Sphere(const vec3& c, float r) : center(c), radius(r) {}

	bool intersect(const Ray& ray, const float& tMin, const float& tMax, float& d) const override {
		vec3 oc = ray.origin - center;
		float a = dot(ray.direction, ray.direction);
		float b = 2.0f * dot(ray.direction, oc);
		float c = dot(oc, oc) - radius * radius;
		float discriminant = b * b - 4 * a * c;
		if (discriminant < 0)
			return false;
		float t1 = (-b - std::sqrt(discriminant)) / (2 * a);
		float t2 = (-b + std::sqrt(discriminant)) / (2 * a);

		if (t1 > tMin && t1 < tMax)
		{
			d = t1;
			return true;
		}
		if (t2 > tMin && t2 < tMax)
		{
			d = t2;
			return true;
		}
		return false;
	}
};

class Plane : public Surface {
public:
	vec3 point, normal;

	Plane(const vec3& p, const vec3& n) : point(p), normal(normalize(n)) {}

	bool intersect(const Ray& ray, const float& tMin, const float& tMax, float& d) const override {
		float denom = dot(normal, ray.direction);
		if (fabs(denom) < 1e-6) return false;

		float t = dot(point - ray.origin, normal) / denom;
		if (t > tMin && t < tMax)
		{
			d = t;
			return true;
		}
		return false;
	}
};

class Scene {
public:
	std::vector<std::shared_ptr<Surface>> objects;

	void addObject(const std::shared_ptr<Surface>& obj) {
		objects.push_back(obj);
	}

	bool trace(const Ray& ray, const float& tMin, const float& tMax) const {
		bool hit = false;
		float closestDistance = INFINITE;

		for (const auto& obj : objects) {
			float distance = 0.0f;
			if (obj->intersect(ray, tMin, tMax, distance)) {
				if (distance < closestDistance) {
					closestDistance = distance;
					hit = true;
				}
			}
		}
		return hit;
	}
};


void render()
{
	//Create our image. We don't want to do this in 
	//the main loop since this may be too slow and we 
	//want a responsive display of our beautiful image.
	//Instead we draw to another buffer and copy this to the 
	//framebuffer using glDrawPixels(...) every refresh
	OutputImage.clear();
	float closest;
	Camera cam(vec3(0, 0, 0), vec3(0, 0, -1), vec3(0, 1, 0),
				-0.1f, 0.1f, -0.1f, 0.1f, 0.1f);
	Scene scene;
	scene.addObject(std::make_shared<Sphere>(vec3(-4.0f, 0.0f, -7.0f), 1.0f));
	scene.addObject(std::make_shared<Sphere>(vec3(0.0f, 0.0f, -7.0f), 2.0f));
	scene.addObject(std::make_shared<Sphere>(vec3(4.0f, 0.0f, -7.0f), 1.0f));
	scene.addObject(std::make_shared<Plane>(vec3(0.0f, -2.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f)));

	for (int j = 0; j < Height; ++j) 
	{
		for (int i = 0; i < Width; ++i) 
		{
			vec3 color = vec3(0.0f, 0.0f, 0.0f); // black color [0,1] in RGB channel

			// ---------------------------------------------------
			// --- Implement your code here to generate the image
			// ---------------------------------------------------

			// in each pixel, check the list of Surface.intersect function
			// if intersect detected, color it white
			// and detect the distance to the intersection from eye(ray)
			Ray ray = cam.generateRay(i, j);
			if (scene.trace(ray, 0.0f, INFINITE))
			{
				color = vec3(1.0f, 1.0f, 1.0f);
			}

			// set the color
			OutputImage.push_back(color.x); // R
			OutputImage.push_back(color.y); // G
			OutputImage.push_back(color.z); // B
		}
	}
}


void resize_callback(GLFWwindow*, int nw, int nh) 
{
	//This is called in response to the window resizing.
	//The new width and height are passed in so we make 
	//any necessary changes:
	Width = nw;
	Height = nh;
	//Tell the viewport to use all of our screen estate
	glViewport(0, 0, nw, nh);

	//This is not necessary, we're just working in 2d so
	//why not let our spaces reflect it?
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(0.0, static_cast<double>(Width)
		, 0.0, static_cast<double>(Height)
		, 1.0, -1.0);

	//Reserve memory for our render so that we don't do 
	//excessive allocations and render the image
	OutputImage.reserve(Width * Height * 3);
	render();
}


int main(int argc, char* argv[])
{
	// -------------------------------------------------
	// Initialize Window
	// -------------------------------------------------

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(Width, Height, "OpenGL Viewer", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	//We have an opengl context now. Everything from here on out 
	//is just managing our window or opengl directly.

	//Tell the opengl state machine we don't want it to make 
	//any assumptions about how pixels are aligned in memory 
	//during transfers between host and device (like glDrawPixels(...) )
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);

	//We call our resize function once to set everything up initially
	//after registering it as a callback with glfw
	glfwSetFramebufferSizeCallback(window, resize_callback);
	resize_callback(NULL, Width, Height);

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		//Clear the screen
		glClear(GL_COLOR_BUFFER_BIT);

		// -------------------------------------------------------------
		//Rendering begins!
		glDrawPixels(Width, Height, GL_RGB, GL_FLOAT, &OutputImage[0]);
		//&OutputImage[0]: viewer window
		//and ends.
		// -------------------------------------------------------------

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();

		//Close when the user hits 'q' or escape
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS
			|| glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
		{
			glfwSetWindowShouldClose(window, GL_TRUE);
		}
	}

	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}

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

	//virtual vec3 getNormal() const = 0;
	virtual vec3 getNormal(const vec3& point) const = 0;

	virtual vec3 getAmbientColor() const = 0;
	virtual vec3 getDiffuseColor() const = 0;
	virtual vec3 getSpecularColor() const = 0;
	virtual float getShininess() const = 0;

	virtual ~Surface() = default;
};

class Sphere : public Surface {
public:
	vec3 center;
	float radius;
	vec3 ka, kd, ks;
	float shininess;

	Sphere(const vec3& c, float r, const vec3& ka, const vec3& kd, const vec3& ks, float shininess)
		: center(c), radius(r), ka(ka), kd(kd), ks(ks), shininess(shininess) {
	}

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

	vec3 getNormal(const vec3& point) const override {
		return normalize(point - center);
	}

	vec3 getAmbientColor() const override { return ka; }
	vec3 getDiffuseColor() const override { return kd; }
	vec3 getSpecularColor() const override { return ks; }
	float getShininess() const override { return shininess; }
};

class Plane : public Surface {
public:
	vec3 point, normal;
	vec3 ka, kd, ks;
	float shininess;

	Plane(const vec3& p, const vec3& n, const vec3& ka, const vec3& kd, const vec3& ks, float shininess)
		: point(p), normal(normalize(n)), ka(ka), kd(kd), ks(ks), shininess(shininess) {
	}

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

	vec3 getNormal(const vec3& point) const override {
		return normal;
	}

	vec3 getAmbientColor() const override { return ka; }
	vec3 getDiffuseColor() const override { return kd; }
	vec3 getSpecularColor() const override { return ks; }
	float getShininess() const override { return shininess; }
};

class PointLight {
public:
	vec3 position;
	vec3 color;

	PointLight(const vec3& position, const vec3& color) : position(position), color(color) {}
};

class Scene {
public:
	std::vector<std::shared_ptr<Surface>> objects;
	PointLight light;

	Scene(const PointLight& light) : light(light) {}

	void addObject(const std::shared_ptr<Surface>& obj) {
		objects.push_back(obj);
	}

	vec3 calculatePhongShading(const Ray& ray, const vec3& hitPoint, const vec3& normal, const Surface* hitObject) const {
		vec3 ambientColor = hitObject->getAmbientColor();
		vec3 diffuseColor = hitObject->getDiffuseColor();
		vec3 specularColor = hitObject->getSpecularColor();
		float shininess = hitObject->getShininess();

		// Ambient Lighting
		vec3 ambient = ambientColor;

		// Shadow Ray
		vec3 lightDir = normalize(light.position - hitPoint);
		Ray shadowRay(hitPoint + lightDir * 1e-4f, lightDir); // Little offset for not reference itself

		bool inShadow = false;
		for (const auto& obj : objects) {
			float shadowT;
			if (obj->intersect(shadowRay, 1e-4f, INFINITY, shadowT)) {
				inShadow = true;
				break;
			}
		}

		if (inShadow) {
			return ambient;  // If in shadow, Duffuse and Specular doesn't occur
		}


		// Diffuse Lighting
		//vec3 lightDir = normalize(light.position - hitPoint);
		float diff = std::max(dot(normal, lightDir), 0.0f);
		vec3 diffuse = diffuseColor * diff;

		// Specular Lighting in Blinn-Phong Model
		vec3 viewDir = normalize(ray.origin - hitPoint);
		vec3 halfVector = normalize(viewDir + lightDir);
		float spec = std::pow(std::max(dot(normal, halfVector), 0.0f), shininess);
		vec3 specular = specularColor * spec;

		// Add All Lightings Up
		vec3 finalColor = ambient + diffuse + specular;
		return finalColor;
	}

	vec3 trace(const Ray& ray, const float& tMin, const float& tMax) const {
		bool hit = false;
		vec3 returnVec = { 0.0f, 0.0f, 0.0f };
		float closestDistance = INFINITE;
		Surface* hitObject = nullptr;

		for (const auto& obj : objects) {
			float distance = 0.0f;
			if (obj->intersect(ray, tMin, tMax, distance)) {
				if (distance < closestDistance) {
					closestDistance = distance;
					hitObject = obj.get();
					hit = true;
				}
			}
		}

		if (hit) {
			vec3 hitPoint = ray.origin + ray.direction * closestDistance;
			vec3 normal = hitObject->getNormal(hitPoint);
			return calculatePhongShading(ray, hitPoint, normal, hitObject);
		}

		return returnVec;
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
	PointLight light({ -4.0, 4.0, -3.0 }, { 0.0, 0.0, 0.0 });
	Scene scene(light);
	scene.addObject(std::make_shared<Plane>(
		vec3(0.0f, -2.0f, 0.0f), //Point
		vec3(0.0f, 1.0f, 0.0f),  //Normal
		vec3(0.2f, 0.2f, 0.2f),  //k ambient
		vec3(1.0f, 1.0f, 1.0f),  //k diffuse
		vec3(0.0f, 0.0f, 0.0f),  //k specular
		0.0f));                  //specular power
	scene.addObject(std::make_shared<Sphere>(
		vec3(-4.0f, 0.0f, -7.0f),//Center
		1.0f, 					 //Radius
		vec3(0.2f, 0.0f, 0.0f),  //k ambient
		vec3(1.0f, 0.0f, 0.0f),  //k diffuse
		vec3(0.0f, 0.0f, 0.0f),  //k specular
		0.0f));					 //specular power
	scene.addObject(std::make_shared<Sphere>(
		vec3(0.0f, 0.0f, -7.0f), //Center
		2.0f,					 //Radius
		vec3(0.0f, 0.2f, 0.0f),  //k ambient
		vec3(0.0f, 0.5f, 0.0f),  //k diffuse
		vec3(0.5f, 0.5f, 0.5f),  //k specular
		32.0f));				 //specular power
	scene.addObject(std::make_shared<Sphere>(
		vec3(4.0f, 0.0f, -7.0f), //Center
		1.0f,					 //Radius
		vec3(0.0f, 0.0f, 0.2f),	 //k ambient
		vec3(0.0f, 0.0f, 1.0f),	 //k diffuse
		vec3(0.0f, 0.0f, 0.0f),	 //k specular
		0.0f));					 //specular power

	const int N = 64; //64 samples for anti-aliasing

	for (int j = 0; j < Height; ++j)
	{
		for (int i = 0; i < Width; ++i)
		{
			vec3 color = vec3(0.0f, 0.0f, 0.0f); // black color [0,1] in RGB channel

			// in each pixel, check the list of Surface.intersect function
			// if intersect detected, color it as trace func returns
			// and detect the distance to the intersection from eye(ray)

			for (int k = 0; k < N; ++k)
			{
				float offsetX = (rand() % 1000) / 1000.0f; // Random Value between [-0.5, 0.5]
				float offsetY = (rand() % 1000) / 1000.0f; // Random Value between [-0.5, 0.5]

				Ray ray = cam.generateRay(i + offsetX, j + offsetY);
				vec3 sampleColor = scene.trace(ray, 0.0f, INFINITE);
				color += sampleColor;
			}
			color /= float(N);

			// set the color, Gamma Correction Applied
			float gamma = 2.2f;
			OutputImage.push_back(pow(color.x, 1.0f / gamma)); // R
			OutputImage.push_back(pow(color.y, 1.0f / gamma)); // G
			OutputImage.push_back(pow(color.z, 1.0f / gamma)); // B
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
	window = glfwCreateWindow(Width, Height, "Q3", NULL, NULL);
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
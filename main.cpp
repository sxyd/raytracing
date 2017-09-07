// Compile using c++ -o main main.cpp -std=c++11 -O3
// ./main

// Shawn Dai

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <stdlib.h>

#include <fstream>
#include <string>
#include <sstream>
using namespace std;

//define class Vec_3, used in ray direction
template<typename T>
class Vec2
{
public:
    Vec2() : x(0), y(0) {}
    Vec2(T xx) : x(xx), y(xx) {}
    Vec2(T xx, T yy) : x(xx), y(yy) {}
    Vec2 operator + (const Vec2 &v) const
    { return Vec2(x + v.x, y + v.y); }
    Vec2 operator / (const T &r) const
    { return Vec2(x / r, y / r); }
    Vec2 operator * (const T &r) const
    { return Vec2(x * r, y * r); }
    Vec2& operator /= (const T &r)
    { x /= r, y /= r; return *this; }
    Vec2& operator *= (const T &r)
    { x *= r, y *= r; return *this; }
    friend std::ostream& operator << (std::ostream &s, const Vec2<T> &v)
    {
        return s << '[' << v.x << ' ' << v.y << ']';
    }
    friend Vec2 operator * (const T &r, const Vec2<T> &v)
    { return Vec2(v.x * r, v.y * r); }
    T x, y;
};

typedef Vec2<float> Vec2f;

template<typename T>
class Vec_3
{
public:
    T x, y, z;
    Vec_3(): x(T(0)), y(T(0)), z(T(0)) {} // p237
    Vec_3(T xx): x(xx), y(xx), z(xx) {}
    Vec_3(T xx, T yy, T zz): x(xx), y(yy), z(zz){}
    Vec_3<T> operator * (const T &f) const
    { return Vec_3<T>(x * f, y * f, z * f);}
    Vec_3<T> operator * (const Vec_3<T> &v) const
    { return Vec_3<T>(x * v.x, y * v.y, z * v.z);}
    T dot(const Vec_3<T> &v) const
    { return x * v.x + y * v.y + z * v.z;}
    Vec_3<T> operator - (const Vec_3<T> &v) const
    { return Vec_3<T>( x - v.x, y - v.y, z - v.z);}
    Vec_3<T> operator + (const  Vec_3<T> &v) const
    { return Vec_3<T>( x + v.x, y + v.y, z + v.z);}
    Vec_3<T>& operator += (const Vec_3<T> &v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    Vec_3<T>& operator *= (const Vec_3<T> &v)
    {
        x *= v.x;
        y *= v.y;
        z *= v.z;
        return *this;
    }
    Vec_3<T> operator - () const
    {
        return Vec_3<T>(-x, -y, -z);
    }
    T length2() const
    {
        return x * x + y * y + z * z;
    }
    T length() const
    {
        return sqrt(length2());
    }
    Vec_3& normal()
    {
        T nor2= length2();
        if (nor2 > 0)
        {
            T nor2_inv= 1 / sqrt(nor2);
            x *= nor2_inv;
            y *= nor2_inv;
            z *= nor2_inv;
        }
        return *this;
    }
    
    Vec_3 crossProduct(const Vec_3<T> &v) const
    { return Vec_3<T>(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }
    
    friend std::ostream & operator << (std::ostream &os, const Vec_3<T> &v)
    {
        os<< "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
};
typedef Vec_3<float> Vec_3f;


inline
float deg2rad(const float &deg)
{ return deg * M_PI / 180; }

inline float modulo(const float &f)
{
    return f - std::floor(f);
}

class Object
{
public:
    virtual ~Object() {}
    // Method to compute the intersection of the object with a ray
    // Returns true if an intersection was found, false otherwise
    // See method implementation in children class for details
    virtual bool intersect(const Vec_3f &, const Vec_3f &, float &, Vec2f &) const = 0;

    virtual Vec_3f getNhit(const Vec_3f &) const = 0;
    virtual float getRlection() const = 0;
    virtual float getTransparency() const = 0;
    virtual Vec_3f getEmissionColor() const = 0;
    virtual Vec_3f getSurfaceColor() const = 0;
    virtual void getSurfaceProperties(const Vec_3f &, const Vec_3f &, const uint32_t &, const Vec2f &, Vec_3f &, Vec2f &) const = 0;
    
    virtual std::string getType() const = 0;
    virtual Vec_3f getCenter() const = 0;
    virtual int isTexture() const = 0;
    
    // Method to compute the surface data such as normal and texture coordnates at the intersection point.
    // See method implementation in children class for details
    std::string type;
    Vec_3f center;
    Vec_3f n;
};


//Define Sphere Class
class Sphere : public Object
{
public:
    Vec_3f center;
    float radius, radius2;
    Vec_3f surfaceColor, emissionColor;
    float transparency, reflection;
    int isTex;
    std::string type = "Sphere";
    Sphere() {};
    Sphere(
           const Vec_3f &c,
           const float &r,
           const Vec_3f &sc,
           const float &refl = 0,
           const float &transp = 0,
           const Vec_3f &ec = 0,
           const int &tex = 0):
    center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
    transparency(transp), reflection(refl), isTex(tex)
    {}
    //Use geometric solution to solve a ray-sphere intersection
    bool intersect(const Vec_3f &rayorigin, const Vec_3f & raydirection, float &t, Vec2f &uv) const
    {
        float t0, t1;
        Vec_3f l = center - rayorigin;
        //Determine whether reverse direction
        float tca = l.dot(raydirection);
        if  (tca < 0) return false;
        //a^2=b^2+c^2
        float dist = l.dot(l) - tca * tca;
        if (dist > radius2) return false;
        float thc = sqrt(radius2 - dist);
        //t0: first intersection distance, t1: second intersection distance
        t0 = tca - thc;
        t1 = tca + thc;
        
        if (t0 > t1) {
            std::swap(t0, t1);
        }
        if (t0 < 0) {
            t0 = t1;
            if (t0 < 0) {
                return false;
            }
        }
        t = t0;
        
        return true;
    }
    
    void getSurfaceProperties(
                            const Vec_3f &hitPoint,
                            const Vec_3f &viewDirection,
                            const uint32_t &triIndex,
                            const Vec2f &uv,
                            Vec_3f &hitNormal,
                            Vec2f &hitTextureCoordinates) const
    {
        hitNormal = hitPoint - center;
        hitNormal.normal();
        // In this particular case, the normal is simular to a point on a unit sphere
        // centred around the origin. We can thus use the normal coordinates to compute
        // the spherical coordinates of Phit.
        // atan2 returns a value in the range [-pi, pi] and we need to remap it to range [0, 1]
        // acosf returns a value in the range [0, pi] and we also need to remap it to the range [0, 1]
        hitTextureCoordinates.x = (1 + atan2(hitNormal.z, hitNormal.x) / M_PI) * 0.5;
        hitTextureCoordinates.y = acosf(hitNormal.y) / M_PI;
    }
    
    Vec_3f getNhit(const Vec_3f &phit) const
    {
        Vec_3f nhit;
        nhit = phit - center;
        return nhit;
    }
    
    float getRlection() const{
        return reflection;
    }
    float getTransparency() const{
        return transparency;
    }
    Vec_3f getEmissionColor() const{
        return emissionColor;
    }
    Vec_3f getSurfaceColor() const{
        return surfaceColor;
    }
    
    std::string getType() const {
        return type;
    }
    
    Vec_3f getCenter() const {
        return center;
    }
    
    int isTexture() const {
        return isTex;
    }
};

//Define Triangle Class
class Triangle : public Object
{
public:
    Vec_3f P0, P1, P2;
    Vec_3f surfaceColor, emissionColor;
    float transparency, reflection;
    std::string type = "Triangle";
    int isTex;
    Triangle() {}
    Triangle(
             const Vec_3f &P0,
             const Vec_3f &P1,
             const Vec_3f &P2,
             const Vec_3f &sc,
             const float &refl = 0,
             const float &transp = 0,
             const Vec_3f &ec = 0,
             const int isTex = 0):
    P0(P0), P1(P1), P2(P2),surfaceColor(sc), reflection(refl), transparency(transp), isTex(isTex)
    {}
    
    float det(Vec_3f a, Vec_3f b, Vec_3f c) const
    {
        return (a.x*b.y*c.z + b.x*c.y*a.z + a.y*b.z*c.x - c.x*b.y*a.z - b.x*a.y*c.z - c.y*b.z*a.x);
    }
    
    bool intersect(const Vec_3f &rayorigin, const Vec_3f & raydirection, float &t, Vec2f &uv) const
    {
        Vec_3f S, E1, E2;
        Vec_3f result;
        Vec_3f ri;
        
        E1 = P0 - P1;
        E2 = P0 - P2;
        S = P0 - rayorigin;
        
        ri.x = det(S, E1, E2);
        ri.y = det(raydirection, S, E2);
        ri.z = det(raydirection, E1, S);
        result = ri * (1 / det(raydirection, E1, E2));
        if (!(result.x > 0 && result.y >= 0 && result.y <= 1 && result.z >=0 &&result.z <= 1 && (result.y + result.z <= 1))) {
            return false;
        }
        t = result.x;
        
        uv.x = result.y;
        uv.y = result.z;
        
        return true;
    }
    
    void getSurfaceProperties(
                              const Vec_3f &hitPoint,
                              const Vec_3f &viewDirection,
                              const uint32_t &triIndex,
                              const Vec2f &uv,
                              Vec_3f &hitNormal,
                              Vec2f &hitTextureCoordinates) const
    {
        Vec2f st0, st1, st2;
        // three situations
        if (hitNormal.x != 0){
            st0.x = P0.z;
            st0.y = P0.y;
            
            st1.x = P1.z;
            st1.y = P1.y;
            
            st2.x = P2.z;
            st2.y = P2.y;
        } else if (hitNormal.y != 0) {
            st0.x = P0.x;
            st0.y = P0.z;
            
            st1.x = P1.x;
            st1.y = P1.z;
            
            st2.x = P2.x;
            st2.y = P2.z;
        } else if (hitNormal.z != 0) {
            st0.x = P0.x;
            st0.y = P0.y;
            
            st1.x = P1.x;
            st1.y = P1.y;
            
            st2.x = P2.x;
            st2.y = P2.y;
        }
        hitTextureCoordinates = (1 - uv.x - uv.y) * st0 + uv.x * st1 + uv.y * st2;
    }
    
    Vec_3f getNhit(const Vec_3f &phit) const
    {
        Vec_3f E1, E2;
        Vec_3f Nhit;
        
        E1 = P0 - P1;
        E2 = P0 - P2;
        
        Nhit.x = E1.y*E2.z - E2.y*E1.z;
        Nhit.y = E1.z*E2.x - E2.z*E1.x;
        Nhit.z = E1.x*E2.y - E2.x*E1.y;
        
        return Nhit;
    }
    float getRlection() const{
        return reflection;
    }
    float getTransparency() const{
        return transparency;
    }
    Vec_3f getEmissionColor() const{
        return Vec_3f(0);
    }
    Vec_3f getSurfaceColor() const{
        return surfaceColor;
    }
    std::string getType() const {
        return type;
    }
    Vec_3f getCenter() const {
        return P0;
    }
    int isTexture() const {
        return isTex;
    }
};

struct IsectInfo
{
    const Object *hitObject = nullptr;
    Vec2f uv;
    uint32_t index = 0;
};

//Define the maximum recursion depth
#define MAX_DEPTH 5

//Calculate the mix value for reflection and refraction
float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}

//Ray Tracing Function: takes a ray (defined by its origin and direction) as argument.
//Through the function, we can know if the ray intersects any of the geometry in the scene.
//If the ray intersects an object, calculate the intersection point and its normal, then shade the point.
//Shading depends on the surface (transparent, reflective, diffuse)
//If the ray intersects an object, then return the color of the object at the intersection point, otherwise return the backgroud color.

Vec_3f trace(
             const Vec_3f &rayorigin,
             const Vec_3f &raydirection,
             const std::vector<std::unique_ptr<Object>> &objects,
             const int &depth
             )
{
    float tnear= INFINITY;
    const Object* object = NULL;
    std::string object_type;
    
    //calculate intersection of this ray with the sphere in the scene
    //find the nearest sphere
    std::vector<std::unique_ptr<Object>>::const_iterator iter = objects.begin();
    
    Vec2f uv;
    for(; iter != objects.end(); ++iter)
    {
        float t = INFINITY;
        
        if((*iter)->intersect(rayorigin, raydirection, t, uv))
        {
            //If the eye point in the sphere
            if(t < tnear)
            {
                tnear = t;
                object = &*(*iter);
            }
        }
    }
    //If there is no intersection, then return backgroud color
    if(!object) return Vec_3f(2); // background color
    
    //Color of ray
    Vec_3f surfaceColor = 0;
    //point of intersect
    Vec_3f phit = rayorigin + raydirection * tnear;
    
    // normal (法线)
    Vec_3f nhit;
    nhit = object->getNhit(phit);
    
    //normalize the normal direction
    nhit.normal();
    //If the normal and the view direction's dot is positive, means the view point inside sphere
    float bias = 1e-3;
    bool inside = false;
    if(raydirection.dot(nhit) > 0)
    {
        nhit = -nhit;
        inside = true;
    }
    
    // add texture
    Vec_3f hitNormal;
    Vec_3f viewDirection;
    uint32_t triIndex;
    Vec2f hitTexCoordinates;
    if (object->getType() == "Sphere") {
        object->getSurfaceProperties(phit, viewDirection, triIndex, uv, hitNormal, hitTexCoordinates);
    } else if (object->getType() == "Triangle") {
        object->getSurfaceProperties(phit, viewDirection, triIndex, uv, nhit, hitTexCoordinates);
    }
    
    //Tackle with relection and refraction
    if((object->getTransparency() > 0 || object->getRlection() > 0) && depth < MAX_DEPTH)
    {
        //Compute fresnel effect
        float facingratio = - raydirection.dot(nhit);
        float fresneleffect = mix (pow(1 - facingratio, 3), 1, 0.1);
        //Compute reflection direction
        Vec_3f reflect_direction = raydirection - nhit * 2 * raydirection.dot(nhit);
        
        reflect_direction.normal();
        Vec_3f next_reflection = trace(phit + nhit * bias, reflect_direction, objects, depth + 1);
        Vec_3f next_refraction = 0;
        //Only if the sphere is transparent, then compute refraction ray
        if(object->getTransparency())
        {
            //judge whether we are inside or outside? ior is the index of two materials
            float ior = 1.1, eta = (inside) ? ior : 1 / ior;
            float cosi = -nhit.dot(raydirection);
            float k = 1 - eta * eta * (1 - cosi * cosi);
            Vec_3f refraction_direction = raydirection * eta + nhit * (eta * cosi - sqrt(k));
            refraction_direction.normal();
            next_refraction = trace(phit - nhit * bias, refraction_direction, objects, depth+1);
        }
        //The surface is a mix of reflection and refraction (if the sphere is transparent)
        surfaceColor = (next_reflection * fresneleffect + next_refraction * (1 - fresneleffect) * object->getTransparency()) * object->getSurfaceColor() ;
    }
    else
    {
        std::vector<std::unique_ptr<Object>>::const_iterator iter = objects.begin();
        for(; iter != objects.end(); ++iter)
        {
            //This is a light
            if((*iter)->getEmissionColor().x > 0)
            {
                Vec_3f transmission = 1;
                Vec_3f lightDirection = (*iter)->getCenter() - phit;
                
                lightDirection.normal();
                //Check whether have an obstacle between light and object, add shadow
                std::vector<std::unique_ptr<Object>>::const_iterator iter2 = objects.begin();
                for(; iter2 != objects.end(); ++iter2)
                {
                    if((*iter2)->getTransparency() >= 0.7)
                    {
                        continue;
                    }
                    if((*iter) != (*iter2))
                    {
                        float t;
                        if((*iter2)->intersect(phit + nhit * bias, lightDirection, t, uv))
                        {
                            transmission = 0;
                            break;
                        }
                    }
                }
                //If nhit and lightDirection's dot is less than 0, then no light.
                surfaceColor += object->getSurfaceColor() * transmission * std::max(float(0), nhit.dot(lightDirection)) * (*iter)->getEmissionColor();

                // compute the pattern
                float angle = deg2rad(0);
                float s = hitTexCoordinates.x * cos(angle) - hitTexCoordinates.y * sin(angle);
                float t = hitTexCoordinates.y * cos(angle) + hitTexCoordinates.x * sin(angle);
                
                float scaleS = 0.3, scaleT = 0.3; // 10 10 0.3 0.3
                if(object->getType() == "Sphere"){
                    scaleS = 10;
                    scaleT = 10;
                }
                float pattern = (modulo(s * scaleS) < 0.5) ^ (modulo(t * scaleT) < 0.5);
                if (object->isTexture() == 1)
                {
                    surfaceColor += transmission * pattern * 0.7 *  std::max(0.f, nhit.dot(lightDirection)) ;
                    
                }
            }
        }
    }
    return surfaceColor + object->getEmissionColor();
}

//Render function, compute each pixel of the image.
void render(const std::vector<std::unique_ptr<Object>> &objects, int j)
{
    unsigned width = 640, height = 480;
    Vec_3f *img = new Vec_3f[width * height], *pixel = img;  //p407
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 60; // angle from top to bottom
    float aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.);
    // Trace all ray
    for(unsigned y = 0; y < height; y++)
    {
        for(unsigned x = 0; x < width; x++, pixel++)
        {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            double m = deg2rad(0);
            double n = deg2rad(0);
            
            double newxx = xx;
            double newyy = yy*cos(m)+(-1)*sin(m);
            double newzz = yy*(-sin(m))+(-1*cos(m));
            
            Vec_3f raydir(newxx*cos(n)+newzz*(-sin(n)), newyy, newxx*sin(n)+newzz*cos(n));
            raydir.normal();
            
            *pixel = trace(Vec_3f(0,0, 30), raydir, objects, 0);
        }
    }
    
    //Save the result
    ostringstream s1;
    int ii = j;
    s1 << ii;
    string s2 = s1.str();
    string str;
    str  = "./" + s2 + ".ppm";
    std::ofstream ofs(str, std::ios::out | std::ios::binary);
    
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for(unsigned i = 0; i < width * height; i++)
    {
        //0,255
        ofs << (unsigned char)(std::min(float(1), img[i].x) * 255) <<
        (unsigned char)(std::min(float(1), img[i].y) * 255) <<
        (unsigned char)(std::min(float(1), img[i].z) * 255);
    }
    ofs.close();
    delete [] img;
}

// read .obj file
void read (vector<Triangle> &triVector1, string filename, int isRflect) {
    string line;

    ifstream myfile (filename); //shape_tetrahedron,
    if (myfile.is_open())
    {
        int raw = 1;
        vector<Vec_3f> pointVector;
        vector<Triangle> triVector;
        
        while ( getline (myfile,line) )
        {
            Vec_3f point;
            string arr[4];
            int ii = 0;
            if(raw == 1)
            {
                raw = 0;
                continue;
            }
            stringstream ssin(line);
            while (ssin.good() && ii < 4){
                ssin >> arr[ii];
                ++ii;
            }
            if (arr[0] == "v")
            {
                point.x = std::stof(arr[1]);
                point.y = std::stof(arr[2]);
                point.z = std::stof(arr[3]);
                pointVector.push_back(point);
            }
            if (arr[0] == "f") {
                int point1 = std::stoi(arr[1]);
                int point2 = std::stoi(arr[2]);
                int point3 = std::stoi(arr[3]);
                
                if(isRflect == 0){
                    triVector.push_back(Triangle(pointVector[point1-1], pointVector[point2-1], pointVector[point3-1], Vec_3f(0.50, 0.00, 0.00), 0, 0));
                } else if(isRflect == 1){
                    triVector.push_back(Triangle(pointVector[point1-1], pointVector[point2-1], pointVector[point3-1], Vec_3f(1, 0.00, 0.00), 1, 0));
                }
            }
        }
        myfile.close();
        triVector1 = triVector;
    }
    else cout << "Unable to open file";
}

//Create a sign including 5 spheres and 1 light (which is also a sphere), then render it.
int main()
{
    std::vector<std::unique_ptr<Object>> objects;
    //argument: position, radius, surfaceColor, reflectivity, transparency, emissionColor
    objects.push_back(std::unique_ptr<Object>(new Sphere(Vec_3f( 15, -8, -5), 2, Vec_3f(0, 3, 0), 1, 0)));// diffuse
    objects.push_back(std::unique_ptr<Object>(new Sphere(Vec_3f( -10, -6, -10), 4, Vec_3f(5), 1, 0)));// diffuse
    objects.push_back(std::unique_ptr<Object>(new Sphere(Vec_3f( -15, -8, -5), 2, Vec_3f(1.00, 0.00, 0.00), 0, 0)));
    objects.push_back(std::unique_ptr<Object>(new Sphere(Vec_3f( 0, -6, -12), 4, Vec_3f(1), 1, 0.8)));
    objects.push_back(std::unique_ptr<Object>(new Sphere(Vec_3f( 10, -6, -10), 4, Vec_3f(0, 0, 1), 0, 0, 0, 1)));
    objects.push_back(std::unique_ptr<Object>(new Sphere(Vec_3f( 0, -9, 5), 1, Vec_3f(0.6, 0.2, 0.6), 0, 0, 0, 1)));
 
    //Triangle,Plane, 2 triangles => squre
    //  down
    objects.push_back(std::unique_ptr<Object>(new Triangle(Vec_3f(-20, -10, -20), Vec_3f(-20, -10, 20), Vec_3f(20, -10, -20), Vec_3f(0.2), 0, 0, 0, 1)));
    objects.push_back(std::unique_ptr<Object>(new Triangle(Vec_3f(20, -10, 20), Vec_3f(20, -10, -20), Vec_3f(-20, -10, 20), Vec_3f(0.2), 0, 0, 0, 1)));
    //  up
    objects.push_back(std::unique_ptr<Object>(new Triangle(Vec_3f(-20, 20, -20), Vec_3f(-20, 20, 20), Vec_3f(20, 20, -20), Vec_3f(0.6), 0, 0)));
    objects.push_back(std::unique_ptr<Object>(new Triangle(Vec_3f(20, 20, 20), Vec_3f(20, 20, -20), Vec_3f(-20, 20, 20), Vec_3f(0.6), 0, 0)));
    //  front
    objects.push_back(std::unique_ptr<Object>(new Triangle(Vec_3f(-20, -10, -20), Vec_3f(-20, 20, -20), Vec_3f(20, 20, -20), Vec_3f(0.2), 0, 0, 0, 1)));
    objects.push_back(std::unique_ptr<Object>(new Triangle(Vec_3f(-20, -10, -20), Vec_3f(20, -10, -20), Vec_3f(20, 20, -20), Vec_3f(0.2), 0, 0, 0, 1)));
    //  left and right
    objects.push_back(std::unique_ptr<Object>(new Triangle(Vec_3f(-20, -10, -20), Vec_3f(-20, -10, 20), Vec_3f(-20, 20, 20), Vec_3f(3), 1, 0)));
    objects.push_back(std::unique_ptr<Object>(new Triangle(Vec_3f(-20, -10, -20), Vec_3f(-20, 20, -20), Vec_3f(-20, 20, 20), Vec_3f(3), 1, 0)));

    objects.push_back(std::unique_ptr<Object>(new Triangle(Vec_3f(20, -10, -20), Vec_3f(20, -10, 20), Vec_3f(20, 50, 20), Vec_3f(4), 1, 0)));
    objects.push_back(std::unique_ptr<Object>(new Triangle(Vec_3f(20, -10, -20), Vec_3f(20, 50, -20), Vec_3f(20, 50, 20), Vec_3f(4), 1, 0)));
    
    vector<Triangle> triVector1;
    read(triVector1, "./cube.obj", 1);
    for(Triangle tri : triVector1)
    {
        objects.push_back(std::unique_ptr<Object>(new Triangle(tri)));
    }
    
    vector<Triangle> triVector2;
    read(triVector2, "./shape_tetrahedron.obj", 0);
    for(Triangle tri2 : triVector2)
    {
        objects.push_back(std::unique_ptr<Object>(new Triangle(tri2)));
    }
    
    int j = 1;
    for (int i = -18; i <= -18; i += 3){
        // Light
        objects.push_back(std::unique_ptr<Object>(new Sphere(Vec_3f(10, 10, 20), 0.5, Vec_3f(0.0, 0.0, 0.0), 0, 0.0, Vec_3f(3))));
        render(objects, j);
        objects.erase(objects.end()-1);
        ++j;
    }
    return 0;
    
}

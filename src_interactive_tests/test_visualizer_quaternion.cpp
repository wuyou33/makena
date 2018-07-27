#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <glm/gtx/norm.hpp>
#include <AntTweakBar.h>

#include "primitives.hpp"
#include "quaternion.hpp"

/** @file src_interactive_tests/test_visualizer_manifold.cpp
 *
 *  @brief visual test rig for Makena::Quaternion
 *
 *  Usage: test_visualizer_manifold
 *
 *  @dependencies
 *    OpenGL 3.3 or later
 *    GLFW3
 *    GLEW
 *    GLM
 *    AntTweakbar 1.16 custom (updated for GLFW3 and Macbook Retina)
 */


static GLFWwindow* window;

static int windowWidth, windowHeight;

static GLuint loadShaders(const char * v_path,const char * f_path);

static void CB_window_size(GLFWwindow *window, int w, int h)
{
    int widthInPixel, heightInPixel;
    glfwGetFramebufferSize(window, &widthInPixel, &heightInPixel);
    glViewport( 0, 0, (GLsizei)widthInPixel, (GLsizei)heightInPixel);
    TwWindowSize(widthInPixel, heightInPixel);
}

static void TwEventMouseButtonGLFW3(
                         GLFWwindow* window, int button, int action, int mods) 
{
    TwEventMouseButtonGLFW(button, action);
}

static void TwEventMousePosGLFW3(GLFWwindow* window, double xpos, double ypos)
{
    TwMouseMotion(int(xpos), int(ypos));
}

static void TwEventMouseWheelGLFW3(
                            GLFWwindow* window, double xoffset, double yoffset)
{
    TwEventMouseWheelGLFW(yoffset);
}

static void TwEventKeyGLFW3(
               GLFWwindow* window, int key, int scancode, int action, int mods)
{
    TwEventKeyGLFW(key, action);
}

static void TwEventCharGLFW3(GLFWwindow* window, int codepoint)
{
    TwEventCharGLFW(codepoint, GLFW_PRESS);
}

static const GLfloat g_vertex_buffer_gcs_axes[] = {

   0.0f, 0.0f, 0.0f,
   1.0f, 0.0f, 0.0f,
   0.0f, 0.0f, 0.0f,
   0.0f, 1.0f, 0.0f,
   0.0f, 0.0f, 0.0f,
   0.0f, 0.0f, 1.0f,

   1.0f, 0.0f, 0.0f,
  10.0f, 0.0f, 0.0f,
   0.0f, 1.0f, 0.0f,
   0.0f,10.0f, 0.0f,
   0.0f, 0.0f, 1.0f,
   0.0f, 0.0f,10.0f,

   0.0f, 0.0f, 0.0f,
 -10.0f, 0.0f, 0.0f,
   0.0f, 0.0f, 0.0f,
   0.0f,-10.0f, 0.0f,
   0.0f, 0.0f, 0.0f,
   0.0f, 0.0f,-10.0f,

};

static const GLfloat g_color_buffer_gcs_axes[] = {

   1.0f, 0.0f, 0.0f,
   1.0f, 0.0f, 0.0f,
   0.0f, 1.0f, 0.0f,
   0.0f, 1.0f, 0.0f,
   0.0f, 0.0f, 1.0f,
   0.0f, 0.0f, 1.0f,

   0.2f, 0.0f, 0.0f,
   0.2f, 0.0f, 0.0f,
   0.0f, 0.2f, 0.0f,
   0.0f, 0.2f, 0.0f,
   0.0f, 0.0f, 0.2f,
   0.0f, 0.0f, 0.2f,

   0.2f, 0.0f, 0.0f,
   0.2f, 0.0f, 0.0f,
   0.0f, 0.2f, 0.0f,
   0.0f, 0.2f, 0.0f,
   0.0f, 0.0f, 0.2f,
   0.0f, 0.0f, 0.2f

};


class GLfloatArrayDelete {public:void operator()(void* x){free(x);}};

using GLfloatArray = std::unique_ptr<GLfloat,GLfloatArrayDelete>;

GLfloatArray makeGLbufferArray(std::vector<Makena::Vec3>& points)
{
    GLfloatArray ar((GLfloat*)malloc(points.size()*sizeof(GLfloat)*3));
    size_t index =0;
    GLfloat* ap = (GLfloat*)ar.get();
    for (auto& p : points) {
        ap[index++] = p.x();
        ap[index++] = p.y();
        ap[index++] = p.z();
    }
    return ar;
}


GLfloatArray makeGLbufferArray(Makena::Vec3& v, size_t numVertices)
{
    GLfloatArray ar((GLfloat*)malloc(numVertices*sizeof(GLfloat)*3));
    size_t index =0;
    GLfloat* ap = (GLfloat*)ar.get();
    for (size_t i = 0; i < numVertices ; i++) {
        ap[index++] = v.x();
        ap[index++] = v.y();
        ap[index++] = v.z();
    }
    return ar;
}



void generateGLBufferFromVecters(
    std::vector<Makena::Vec3>& vectors,
    GLfloat* buf
) {
    GLfloat* fp = buf;
    for (auto& v : vectors) {
        *(fp++) = v.x();
        *(fp++) = v.y();
        *(fp++) = v.z();
    }
}


void generateGLBufferFromVecters4(
    std::vector<Makena::Vec3>& vectors,
    float w,
    GLfloat* buf
) {
    GLfloat* fp = buf;
    for (auto& v : vectors) {
        *(fp++) = v.x();
        *(fp++) = v.y();
        *(fp++) = v.z();
        *(fp++) = w;
    }
}



void generateVerticesFromQuaternion_GL_TRIANGLES(
    std::vector<glm::quat*>&              quats,
    std::vector<Makena::Vec3>& vertices,
    std::vector<Makena::Vec3>& colors,
    std::vector<Makena::Vec3>& normals
) {
    vertices.clear();
    colors.clear();
    normals.clear();

    Makena::Vec3 origin(0.0, 0.0, 0.0);
    Makena::Vec3 unitX(1.0, 0.0, 0.0);
    Makena::Vec3 unitY(0.0, 1.0, 0.0);
    Makena::Vec3 unitZ(0.0, 0.0, 1.0);

    Makena::Vec3 black(0.0, 0.0, 0.0);
    Makena::Vec3 red(1.0, 0.0, 0.0);
    Makena::Vec3 green(0.0, 1.0, 0.0);
    Makena::Vec3 blue(0.0, 0.0, 1.0);

    Makena::Vec3 normX1(1.0, 0.0, 0.0);
    Makena::Vec3 normX2(-1.0, 0.0, 0.0);
    Makena::Vec3 normY1(0.0, 1.0, 0.0);
    Makena::Vec3 normY2(0.0, -1.0, 0.0);
    Makena::Vec3 normZ1(0.0, 0.0, 1.0);
    Makena::Vec3 normZ2(0.0, 0.0, -1.0);

    for (auto& q : quats) {

        Makena::Quaternion tQ(q->w, q->x, q->y, q->z);
//        tQ.normalize();
        auto Mrot    = tQ.rotationMatrix();
        auto runitX  = Mrot * unitX;
        auto runitY  = Mrot * unitY;
        auto runitZ  = Mrot * unitZ;
        auto rnormX1 = Mrot * normX1;
        auto rnormX2 = Mrot * normX2;
        auto rnormY1 = Mrot * normY1;
        auto rnormY2 = Mrot * normY2;
        auto rnormZ1 = Mrot * normZ1;
        auto rnormZ2 = Mrot * normZ2;

        // Face 1 & 2
        vertices.push_back(origin);
        vertices.push_back(runitX);
        vertices.push_back(runitY);

        colors.push_back(black);
        colors.push_back(red);
        colors.push_back(green);

        normals.push_back(normZ1);
        normals.push_back(normZ1);
        normals.push_back(normZ1);

        vertices.push_back(origin);
        vertices.push_back(runitY);
        vertices.push_back(runitX);

        colors.push_back(black);
        colors.push_back(green);
        colors.push_back(red);

        normals.push_back(normZ2);
        normals.push_back(normZ2);
        normals.push_back(normZ2);

        // Face 3 & 4
        vertices.push_back(origin);
        vertices.push_back(runitY);
        vertices.push_back(runitZ);

        colors.push_back(black);
        colors.push_back(green);
        colors.push_back(blue);

        normals.push_back(normX1);
        normals.push_back(normX1);
        normals.push_back(normX1);

        vertices.push_back(origin);
        vertices.push_back(runitZ);
        vertices.push_back(runitY);

        colors.push_back(black);
        colors.push_back(blue);
        colors.push_back(green);

        normals.push_back(normX2);
        normals.push_back(normX2);
        normals.push_back(normX2);

        // Face 5 & 6
        vertices.push_back(origin);
        vertices.push_back(runitZ);
        vertices.push_back(runitX);

        colors.push_back(black);
        colors.push_back(blue);
        colors.push_back(red);

        normals.push_back(normY1);
        normals.push_back(normY1);
        normals.push_back(normY1);

        vertices.push_back(origin);
        vertices.push_back(runitX);
        vertices.push_back(runitZ);

        colors.push_back(black);
        colors.push_back(red);
        colors.push_back(blue);

        normals.push_back(normY2);
        normals.push_back(normY2);
        normals.push_back(normY2);

    }    
}


void generateVerticesFromQuaternion_GL_LINES(
    Makena::Quaternion&        quat,
    std::vector<Makena::Vec3>& vertices,
    std::vector<Makena::Vec3>& colors
) {
    vertices.clear();
    colors.clear();

    Makena::Vec3 origin(0.0, 0.0, 0.0);
    Makena::Vec3 unitX(1.0, 0.0, 0.0);
    Makena::Vec3 unitY(0.0, 1.0, 0.0);
    Makena::Vec3 unitZ(0.0, 0.0, 1.0);
    Makena::Vec3 red(1.0, 0.3, 0.3);
    Makena::Vec3 green(0.3, 1.0, 0.3);
    Makena::Vec3 blue(0.3, 0.3, 1.0);

    auto& Mrot = quat.rotationMatrix();    

    auto rotX = Mrot * unitX;
    auto rotY = Mrot * unitY;
    auto rotZ = Mrot * unitZ;

    vertices.push_back(origin);
    vertices.push_back(rotX);
    colors.push_back(red);
    colors.push_back(red);

    vertices.push_back(origin);
    vertices.push_back(rotY);
    colors.push_back(green);
    colors.push_back(green);

    vertices.push_back(origin);
    vertices.push_back(rotZ);
    colors.push_back(blue);
    colors.push_back(blue);
}


void generateVerticesFromQuaternion_GL_POINTS(
    Makena::Quaternion&        quat,
    std::vector<Makena::Vec3>& vertices,
    std::vector<Makena::Vec3>& colors
) {
    vertices.clear();
    colors.clear();

    Makena::Vec3 unitX(1.0, 0.0, 0.0);
    Makena::Vec3 unitY(0.0, 1.0, 0.0);
    Makena::Vec3 unitZ(0.0, 0.0, 1.0);
    Makena::Vec3 red(1.0, 0.3, 0.3);
    Makena::Vec3 green(0.3, 1.0, 0.3);
    Makena::Vec3 blue(0.3, 0.3, 1.0);

    auto& Mrot = quat.rotationMatrix();    

    auto rotX = Mrot * unitX;
    auto rotY = Mrot * unitY;
    auto rotZ = Mrot * unitZ;

    vertices.push_back(rotX);
    colors.push_back(red);

    vertices.push_back(rotY);
    colors.push_back(green);

    vertices.push_back(rotZ);
    colors.push_back(blue);
}


int main( int argc, char* argv[])
{

    if (argc != 2) {
        std::cerr << "Usage: test_visualizer_quaternion "
                     "[num of quaternions] "
                     "If 2 is given, then this works as an interpolator mode";
        exit(1);
    }

    int numQuats = atoi(argv[1]);

    if (numQuats < 1 || numQuats > 20) {
        std::cerr << "Usage: test_visualizer_quaternion "
                     "[num of quaternions] "
                     "If 2 is given, then this works as an interpolator mode";
        exit(1);
    }        

    // Initialise GLFW
    if (!glfwInit()) {
        std::cerr << "glfwInit() failed.\n";
        return 1;
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window = glfwCreateWindow(
                         1024, 768, "Test Visualiser Quaternion", NULL, NULL);

    if ( window == NULL ) {
        std::cerr << "glfwCreateWindow failed\n";
        glfwTerminate();
        return 1;
    }

    glfwMakeContextCurrent(window);

    int widthInPixel, heightInPixel;
    glfwGetFramebufferSize(window, &widthInPixel, &heightInPixel);

    int width, height;
    glfwGetWindowSize(window, &width, &height);

    glfwSetWindowSizeCallback(window, CB_window_size);
    glfwGetWindowSize(window, &windowWidth, &windowHeight);

    glewExperimental = true;

    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW\n";;
        glfwTerminate();
        return 1;
    }

    TwInit(TW_OPENGL_CORE, NULL);
    TwWindowSize(widthInPixel, heightInPixel);
    TwSetPixelRatio(
         float(widthInPixel)/(float)width, float(heightInPixel)/(float)height);
    TwBar* twBar01 = TwNewBar("Control");

    TwDefine(" GLOBAL help='Control Panel' ");

    double distance = 7.0; 
    TwAddVarRW(twBar01, "distance", TW_TYPE_DOUBLE, &distance, 
        " label='Camera Distance' min=2.0 max=20.0 step=0.1 "
        "keyIncr=s keyDecr=S help='Rotation speed (turns/second)' ");

    glm::quat orientation;
    TwAddVarRW(twBar01, "Model Rotation", TW_TYPE_QUAT4F, &orientation, 
                                                 " showval=true opened=true ");

    TwSetParam(twBar01, NULL, "refresh",  TW_PARAM_CSTRING, 1, "0.1");
    TwSetParam(twBar01, NULL, "position", TW_PARAM_CSTRING, 1, " 10 10 ");

    std::vector<TwBar*>     twBars;
    std::vector<glm::quat*> twQuats;


    for (int i = 0; i < numQuats; i++) {
        char title[10];
        sprintf(title, "Quat %d", i+1);
        TwBar*     tp = TwNewBar(title);
        glm::quat* qp = new glm::quat(1.0, 0.0, 0.0, 0.0);
        TwAddVarRW(tp, title, TW_TYPE_QUAT4F, qp," showval=true opened=true ");
        twBars.push_back (tp);
        twQuats.push_back(qp);
    }

    glfwSetMouseButtonCallback(window,(GLFWmousebuttonfun)
                                                     TwEventMouseButtonGLFW3);
    glfwSetCursorPosCallback  (window,(GLFWcursorposfun)
                                                     TwEventMousePosGLFW3);
    glfwSetScrollCallback     (window,(GLFWscrollfun)TwEventMouseWheelGLFW3);
    glfwSetKeyCallback        (window,(GLFWkeyfun)   TwEventKeyGLFW3);
    glfwSetCharCallback       (window,(GLFWcharfun)  TwEventCharGLFW3);

    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    glfwPollEvents();

    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS); 
    glEnable(GL_CULL_FACE);

    GLuint vertexBufferGCSAxes;
    glGenBuffers(1, &vertexBufferGCSAxes);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBufferGCSAxes);
    glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(g_vertex_buffer_gcs_axes),
        g_vertex_buffer_gcs_axes,
        GL_STATIC_DRAW
    );

    GLuint colorBufferGCSAxes;
    glGenBuffers(1, &colorBufferGCSAxes);
    glBindBuffer(GL_ARRAY_BUFFER, colorBufferGCSAxes);
    glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(g_color_buffer_gcs_axes),
        g_color_buffer_gcs_axes,
        GL_STATIC_DRAW
    );

    Makena::Vec3 dummyVert(1.0, 0.0, 0.0);
    auto vertex_buffer_quaternions = makeGLbufferArray(dummyVert,
        numQuats * 
        6 /* faces */   *
        3 /* vertices */  );
    GLuint vertexBufferQuaternions;
    glGenBuffers(1, &vertexBufferQuaternions);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBufferQuaternions);
    glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(GLfloat)*(numQuats)*6*3*3,
        vertex_buffer_quaternions.get(),
        GL_DYNAMIC_DRAW
    );

    auto color_buffer_quaternions = makeGLbufferArray(dummyVert,
        numQuats * 
        6 /* faces */   *
        4 /* rgba */  );
    GLuint colorBufferQuaternions;
    glGenBuffers(1, &colorBufferQuaternions);
    glBindBuffer(GL_ARRAY_BUFFER, colorBufferQuaternions);
    glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(GLfloat)*(numQuats)*6*3*4,
        color_buffer_quaternions.get(),
        GL_DYNAMIC_DRAW
    );

    auto normal_buffer_quaternions = makeGLbufferArray(dummyVert,
        numQuats * 
        6 /* faces */   *
        3 /* vertices */  );
    GLuint normalBufferQuaternions;
    glGenBuffers(1, &normalBufferQuaternions);
    glBindBuffer(GL_ARRAY_BUFFER, normalBufferQuaternions);
    glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(GLfloat)*(numQuats)*6*3*3,
        normal_buffer_quaternions.get(),
        GL_DYNAMIC_DRAW
    );

    auto vertex_buffer_average = makeGLbufferArray(dummyVert, 3 * 2);
    GLuint vertexBufferAverage;
    glGenBuffers(1, &vertexBufferAverage);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBufferAverage);
    glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(GLfloat)*3*2*3,
        vertex_buffer_average.get(),
        GL_DYNAMIC_DRAW
    );

    auto color_buffer_average = makeGLbufferArray(dummyVert, 3 * 2);
    GLuint colorBufferAverage;
    glGenBuffers(1, &colorBufferAverage);
    glBindBuffer(GL_ARRAY_BUFFER, colorBufferAverage);
    glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(GLfloat)*3*2*3,
        color_buffer_average.get(),
        GL_DYNAMIC_DRAW
    );

    auto vertex_buffer_average_pts = makeGLbufferArray(dummyVert, 3);
    GLuint vertexBufferAveragePts;
    glGenBuffers(1, &vertexBufferAveragePts);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBufferAveragePts);
    glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(GLfloat)*3*3,
        vertex_buffer_average_pts.get(),
        GL_DYNAMIC_DRAW
    );

    auto color_buffer_average_pts = makeGLbufferArray(dummyVert, 3);
    GLuint colorBufferAveragePts;
    glGenBuffers(1, &colorBufferAveragePts);
    glBindBuffer(GL_ARRAY_BUFFER, colorBufferAveragePts);
    glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(GLfloat)*3*3,
        color_buffer_average_pts.get(),
        GL_DYNAMIC_DRAW
    );
 
    GLuint SID_POL = loadShaders(
        "src_interactive_tests/shaders/vert_pol.glsl",
        "src_interactive_tests/shaders/frag_pol.glsl" 
    );

    GLuint SID_PT = loadShaders(
        "src_interactive_tests/shaders/vert_pt.glsl",
        "src_interactive_tests/shaders/frag_pt.glsl" 
    );

    GLuint SID_POL_MVP       = glGetUniformLocation(SID_POL, "MVP");
    GLuint SID_POL_M         = glGetUniformLocation(SID_POL, "M");
    GLuint SID_POL_V         = glGetUniformLocation(SID_POL, "V");
    GLuint SID_POL_lightWCS  = glGetUniformLocation(SID_POL, "lightWCS");

    GLuint SID_POL_vertexLCS = glGetAttribLocation(SID_POL, "vertexLCS");
    GLuint SID_POL_vertexRGBA= glGetAttribLocation(SID_POL, "vertexRGBA");
    GLuint SID_POL_normalLCS = glGetAttribLocation(SID_POL, "normalLCS");

    GLuint SID_PT_MVP       = glGetUniformLocation(SID_PT, "MVP");
    GLuint SID_PT_M         = glGetUniformLocation(SID_PT, "M");
    GLuint SID_PT_lightWCS  = glGetUniformLocation(SID_PT, "lightWCS");

    GLuint SID_PT_vertexLCS = glGetAttribLocation(SID_PT, "vertexLCS");
    GLuint SID_PT_vertexRGB = glGetAttribLocation(SID_PT, "vertexRGB");

    double alpha = 1.0;
    bool   goingUp = false;
    long loopCount = 0;

    do {
        if (goingUp) {
            alpha +=     0.01;
            if (alpha >= 1.0)
                goingUp = false;
        }
        else {
            alpha -=     0.01;
            if (alpha <= 0.0)
                goingUp = true;

        }
        std::vector<Makena::Quaternion> quats;
        std::vector<double>                        weights;
        for (size_t i = 0; i < twQuats.size(); i++) {
            auto qp = twQuats[i];
            Makena::Quaternion q(qp->w, qp->x, qp->y, qp->z);
            q.normalize();
            quats.push_back(q);         
            if (numQuats != 2) {
                weights.push_back(1.0/(numQuats));
            }
            else {
                if (i == 0)
                    weights.push_back(alpha);
                else if (i == 1)
                    weights.push_back((1.0 - alpha));
            }
        }
 
        auto Qavg = Makena::Quaternion::average(quats, weights);
        if (loopCount++ % 100 == 0) {
            for (size_t i = 0; i < twQuats.size(); i++) {
                std::cerr << "Q[" << i+1 << "]: " 
                    << twQuats[i]->w << ","
                    << twQuats[i]->x << ","
                    << twQuats[i]->y << ","
                    << twQuats[i]->z << "\n";
            }
            std::cerr << "AVG: " << Qavg.s() << "," << 
                                    Qavg.x() << "," << 
                                    Qavg.y() << "," << 
                                    Qavg.z() << "\n";
        }

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        std::vector<Makena::Vec3> vecVerticesQuats;
        std::vector<Makena::Vec3> vecColorsQuats;
        std::vector<Makena::Vec3> vecNormalsQuats;

        generateVerticesFromQuaternion_GL_TRIANGLES(
                   twQuats, vecVerticesQuats, vecColorsQuats, vecNormalsQuats);

        generateGLBufferFromVecters(
                            vecVerticesQuats, vertex_buffer_quaternions.get());

        generateGLBufferFromVecters4(
                          vecColorsQuats, 0.8, color_buffer_quaternions.get());

        generateGLBufferFromVecters(
                             vecNormalsQuats, normal_buffer_quaternions.get());

        std::vector<Makena::Vec3> vecVerticesAvg;
        std::vector<Makena::Vec3> vecColorsAvg;

        generateVerticesFromQuaternion_GL_LINES(
                                           Qavg, vecVerticesAvg, vecColorsAvg);

        generateGLBufferFromVecters(
                                  vecVerticesAvg, vertex_buffer_average.get());

        generateGLBufferFromVecters(vecColorsAvg, color_buffer_average.get());


        std::vector<Makena::Vec3> vecVerticesAvgPts;
        std::vector<Makena::Vec3> vecColorsAvgPts;

        generateVerticesFromQuaternion_GL_POINTS(
                                     Qavg, vecVerticesAvgPts, vecColorsAvgPts);

        generateGLBufferFromVecters(
                           vecVerticesAvgPts, vertex_buffer_average_pts.get());

        generateGLBufferFromVecters(
                              vecColorsAvgPts, color_buffer_average_pts.get());



        glm::vec3 scaleFactor(1.0f, 1.0f, 1.0f);
        glm::mat4 Mscale = glm::scale( glm::mat4(), scaleFactor );
                                   
        glm::quat Qrot(0.0f, 1.0f, 0.0f, 0.0f);
	glm::mat4 Mrot   = glm::mat4_cast(orientation);

	glm::mat4 Mtrans = glm::translate(
                                    glm::mat4(), glm::vec3(0.0f, 0.0f, 0.0f));
	glm::mat4 M      = Mtrans * Mrot * Mscale;

	glm::mat4 V      = glm::lookAt(glm::vec3( 0, 0, distance ), // Cam pos
                                       glm::vec3( 0, 0, 0 ),  // and looks here
                                       glm::vec3( 0, 1, 0 ) );// Head is up

	glm::mat4 P      = glm::perspective(
                               glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
 
        glm::mat4 MVP    = P * V * M;

	glm::vec3 lightPos = glm::vec3(4.0, 4.0, 4.0);

        {
            glDisable(GL_BLEND);
            glUseProgram(SID_PT);



            glUniformMatrix4fv(SID_PT_MVP, 1, GL_FALSE, &MVP[0][0] );
            glUniformMatrix4fv(SID_PT_M,   1, GL_FALSE, &M  [0][0] );
            glUniform3f(SID_PT_lightWCS, lightPos.x, lightPos.y, lightPos.z);

            glEnableVertexAttribArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, vertexBufferGCSAxes);
            glVertexAttribPointer(
                         SID_PT_vertexLCS, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, colorBufferGCSAxes);
            glVertexAttribPointer(
                         SID_PT_vertexRGB, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glDrawArrays(GL_LINES, 0, 6*3);

            glDisableVertexAttribArray(1);
            glDisableVertexAttribArray(0);

            glEnableVertexAttribArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, vertexBufferAverage);
            glBufferData(
                GL_ARRAY_BUFFER,
                sizeof(GLfloat)*3*2*3,
                vertex_buffer_average.get(),
                GL_DYNAMIC_DRAW
            );
            glVertexAttribPointer(
                         SID_PT_vertexLCS, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, colorBufferAverage);
            glBufferData(
                GL_ARRAY_BUFFER,
                sizeof(GLfloat)*3*2*3,
                color_buffer_average.get(),
                GL_DYNAMIC_DRAW
            );
            glVertexAttribPointer(
                         SID_PT_vertexRGB, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glDrawArrays(GL_LINES, 0, 6);

            glDisableVertexAttribArray(1);
            glDisableVertexAttribArray(0);


            glEnableVertexAttribArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, vertexBufferAveragePts);
            glBufferData(
                GL_ARRAY_BUFFER,
                sizeof(GLfloat)*3*3,
                vertex_buffer_average_pts.get(),
                GL_DYNAMIC_DRAW
            );
            glVertexAttribPointer(
                         SID_PT_vertexLCS, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, colorBufferAveragePts);
            glBufferData(
                GL_ARRAY_BUFFER,
                sizeof(GLfloat)*3*3,
                color_buffer_average_pts.get(),
                GL_DYNAMIC_DRAW
            );
            glVertexAttribPointer(
                         SID_PT_vertexRGB, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glPointSize(5.0);
            glDrawArrays(GL_POINTS, 0, 3);

            glDisableVertexAttribArray(1);
            glDisableVertexAttribArray(0);

        }


        {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            glUseProgram(SID_POL);

            glUniformMatrix4fv(SID_POL_MVP, 1, GL_FALSE, &MVP[0][0] );
            glUniformMatrix4fv(SID_POL_M,   1, GL_FALSE, &M  [0][0] );
            glUniformMatrix4fv(SID_POL_V,   1, GL_FALSE, &V  [0][0] );
            glUniform3f(SID_POL_lightWCS, lightPos.x, lightPos.y, lightPos.z);

            glEnableVertexAttribArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, vertexBufferQuaternions);
            glBufferData(
                GL_ARRAY_BUFFER,
                sizeof(GLfloat)*numQuats*6*3*3,
                vertex_buffer_quaternions.get(),
                GL_DYNAMIC_DRAW
            );
            glVertexAttribPointer(
                         SID_POL_vertexLCS, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, colorBufferQuaternions);
            glBufferData(
                GL_ARRAY_BUFFER,
                sizeof(GLfloat)*numQuats*6*3*4,
                color_buffer_quaternions.get(),
                GL_DYNAMIC_DRAW
            );
            glVertexAttribPointer(
                         SID_POL_vertexRGBA, 4, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glEnableVertexAttribArray(2);
            glBindBuffer(GL_ARRAY_BUFFER, normalBufferQuaternions);
            glBufferData(
                GL_ARRAY_BUFFER,
                sizeof(GLfloat)*numQuats*6*3*3,
                normal_buffer_quaternions.get(),
                GL_DYNAMIC_DRAW
            );

            glVertexAttribPointer(
                         SID_POL_normalLCS, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glDrawArrays(GL_TRIANGLES, 0, 
                numQuats * 
                6 /* faces */   *
                3 /* vertices */  );


            glDisableVertexAttribArray(2);
            glDisableVertexAttribArray(1);
            glDisableVertexAttribArray(0);
        }

        TwDraw();

        glfwSwapBuffers(window);

        glfwPollEvents();

    } while ( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
              glfwWindowShouldClose(window) == 0                    );

    TwTerminate();

    for (auto tp: twQuats) {
        delete tp;
    }

    glDeleteBuffers(1, &vertexBufferQuaternions);
    glDeleteBuffers(1, &colorBufferQuaternions);
    glDeleteBuffers(1, &normalBufferQuaternions);

    glDeleteBuffers(1, &vertexBufferGCSAxes);
    glDeleteBuffers(1, &colorBufferGCSAxes);
    glDeleteVertexArrays(1, &VertexArrayID);

    glDeleteProgram(SID_PT);

    glfwTerminate();

    return 0;
}


static GLuint loadShaders(const char * v_path,const char * f_path)
{
    GLuint VertexShaderID   = glCreateShader(GL_VERTEX_SHADER);
    GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
    std::string   VertexShaderCode;
    std::ifstream VertexShaderStream(v_path, std::ios::in);

    if (VertexShaderStream.is_open()) {

        std::string Line = "";

        while (getline(VertexShaderStream, Line)) {

            VertexShaderCode += "\n" + Line;
        }

        VertexShaderStream.close();

    }
    else {
        std::cerr << "failed to open [%s]\n";
        return 0;
    }

    std::string   FragmentShaderCode;
    std::ifstream FragmentShaderStream(f_path, std::ios::in);

    if(FragmentShaderStream.is_open()){

        std::string Line = "";

        while (getline(FragmentShaderStream, Line)) {

            FragmentShaderCode += "\n" + Line;

        }

        FragmentShaderStream.close();
    }

    GLint Result = GL_FALSE;
    int   InfoLogLength;

    char const * VertexSourcePointer = VertexShaderCode.c_str();
    glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
    glCompileShader(VertexShaderID);

    glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);

    if ( InfoLogLength > 0 ) {

        std::vector<char> VertexShaderErrorMessage(InfoLogLength+1);
        glGetShaderInfoLog(
            VertexShaderID,
            InfoLogLength,
            NULL,
            &VertexShaderErrorMessage[0]
        );
        std::cerr << "[" << v_path << "]: ";
        std::cerr << &VertexShaderErrorMessage[0] << "\n";
    }

    char const * FragmentSourcePointer = FragmentShaderCode.c_str();
    glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
    glCompileShader(FragmentShaderID);

    glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if ( InfoLogLength > 0 ) {

        std::vector<char> FragmentShaderErrorMessage(InfoLogLength+1);
        glGetShaderInfoLog(
            FragmentShaderID,
            InfoLogLength,
            NULL,
            &FragmentShaderErrorMessage[0]
        );
        std::cerr << "[" << f_path << "]: ";
        std::cerr << &FragmentShaderErrorMessage[0] << "\n";;
    }

    GLuint ProgramID = glCreateProgram();
    glAttachShader(ProgramID, VertexShaderID);
    glAttachShader(ProgramID, FragmentShaderID);
    glLinkProgram(ProgramID);
    glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
    glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);

    if ( InfoLogLength > 0 ) {
        std::vector<char> ProgramErrorMessage(InfoLogLength+1);
        glGetProgramInfoLog(
                      ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
        std::cerr << "Linker: ";
        std::cerr << &ProgramErrorMessage[0] << "\n";

    }


    glDetachShader(ProgramID, VertexShaderID);
    glDetachShader(ProgramID, FragmentShaderID);
  
    glDeleteShader(VertexShaderID);
    glDeleteShader(FragmentShaderID);

    return ProgramID;
}


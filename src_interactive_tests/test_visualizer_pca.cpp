#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <glm/gtx/norm.hpp>
#include <AntTweakBar.h>

#include "primitives.hpp"

/** @file src_interactive_tests/test_visualizer_pca.cpp
 *
 *  @brief visual test rig for Makena::findPrincipalComponents()
 *
 *  Usage: test_visualizer_pca <test_input_file>
 *
 *
 *  Accompanying test pattern generator is util/generate_random_3d_plots.py
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

static void CB_window_size(GLFWwindow *window, int w, int h){

    int widthInPixel, heightInPixel;
    glfwGetFramebufferSize(window, &widthInPixel, &heightInPixel);
    glViewport( 0, 0, (GLsizei)widthInPixel, (GLsizei)heightInPixel);
    TwWindowSize(widthInPixel, heightInPixel);

}

void CB_cursor_pos(GLFWwindow *window, double x, double y){
    std::cerr << "CB cursor pos: " << x << "," << y << "\n";
}

static void TwEventMouseButtonGLFW3(GLFWwindow* window, int button, int action, int mods) {
    TwEventMouseButtonGLFW(button, action);
}

static void TwEventMousePosGLFW3(GLFWwindow* window, double xpos, double ypos) {
    TwMouseMotion(int(xpos), int(ypos));
}

static void TwEventMouseWheelGLFW3(GLFWwindow* window, double xoffset, double yoffset) {
    TwEventMouseWheelGLFW(yoffset);
}

static void TwEventKeyGLFW3(GLFWwindow* window, int key, int scancode, int action, int mods) {
    TwEventKeyGLFW(key, action);
}

static void TwEventCharGLFW3(GLFWwindow* window, int codepoint) {
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


static const GLfloat g_vertex_buffer_data[] = { 
  -1.0f,-1.0f,-1.0f,
  -1.0f,-1.0f, 1.0f,
  -1.0f, 1.0f, 1.0f,

  -1.0f,-1.0f,-1.0f,
  -1.0f, 1.0f, 1.0f,
  -1.0f, 1.0f,-1.0f,

   1.0f, 1.0f, 1.0f,
   1.0f,-1.0f,-1.0f,
   1.0f, 1.0f,-1.0f,

   1.0f,-1.0f,-1.0f,
   1.0f, 1.0f, 1.0f,
   1.0f,-1.0f, 1.0f,

   1.0f,-1.0f, 1.0f,
  -1.0f,-1.0f, 1.0f,
  -1.0f,-1.0f,-1.0f,

   1.0f,-1.0f, 1.0f,
  -1.0f,-1.0f,-1.0f,
   1.0f,-1.0f,-1.0f,

   1.0f, 1.0f, 1.0f,
   1.0f, 1.0f,-1.0f,
  -1.0f, 1.0f,-1.0f,

   1.0f, 1.0f, 1.0f,
  -1.0f, 1.0f,-1.0f,
  -1.0f, 1.0f, 1.0f,

   1.0f, 1.0f,-1.0f,
  -1.0f,-1.0f,-1.0f,
  -1.0f, 1.0f,-1.0f,

   1.0f, 1.0f,-1.0f,
   1.0f,-1.0f,-1.0f,
  -1.0f,-1.0f,-1.0f,

  -1.0f, 1.0f, 1.0f,
  -1.0f,-1.0f, 1.0f,
   1.0f,-1.0f, 1.0f,

   1.0f, 1.0f, 1.0f,
  -1.0f, 1.0f, 1.0f,
   1.0f,-1.0f, 1.0f
};


static const GLfloat g_normal_buffer_data[] = { 
  -1.0f, 0.0f, 0.0f,
  -1.0f, 0.0f, 0.0f,
  -1.0f, 0.0f, 0.0f,

  -1.0f, 0.0f, 0.0f,
  -1.0f, 0.0f, 0.0f,
  -1.0f, 0.0f, 0.0f,

   1.0f, 0.0f, 0.0f,
   1.0f, 0.0f, 0.0f,
   1.0f, 0.0f, 0.0f,

   1.0f, 0.0f, 0.0f,
   1.0f, 0.0f, 0.0f,
   1.0f, 0.0f, 0.0f,

   0.0f,-1.0f, 0.0f,
   0.0f,-1.0f, 0.0f,
   0.0f,-1.0f, 0.0f,

   0.0f,-1.0f, 0.0f,
   0.0f,-1.0f, 0.0f,
   0.0f,-1.0f, 0.0f,

   0.0f, 1.0f, 0.0f,
   0.0f, 1.0f, 0.0f,
   0.0f, 1.0f, 0.0f,

   0.0f, 1.0f, 0.0f,
   0.0f, 1.0f, 0.0f,
   0.0f, 1.0f, 0.0f,

   0.0f, 0.0f,-1.0f,
   0.0f, 0.0f,-1.0f,
   0.0f, 0.0f,-1.0f,

   0.0f, 0.0f,-1.0f,
   0.0f, 0.0f,-1.0f,
   0.0f, 0.0f,-1.0f,

   0.0f, 0.0f, 1.0f,
   0.0f, 0.0f, 1.0f,
   0.0f, 0.0f, 1.0f,

   0.0f, 0.0f, 1.0f,
   0.0f, 0.0f, 1.0f,
   0.0f, 0.0f, 1.0f,
};


static const GLfloat g_color_buffer_data[] = { 
  0.583f,  0.771f,  0.014f,
  0.609f,  0.115f,  0.436f,
  0.327f,  0.483f,  0.844f,
  0.822f,  0.569f,  0.201f,
  0.435f,  0.602f,  0.223f,
  0.310f,  0.747f,  0.185f,
  0.597f,  0.770f,  0.761f,
  0.559f,  0.436f,  0.730f,
  0.359f,  0.583f,  0.152f,
  0.483f,  0.596f,  0.789f,
  0.559f,  0.861f,  0.639f,
  0.195f,  0.548f,  0.859f,
  0.014f,  0.184f,  0.576f,
  0.771f,  0.328f,  0.970f,
  0.406f,  0.615f,  0.116f,
  0.676f,  0.977f,  0.133f,
  0.971f,  0.572f,  0.833f,
  0.140f,  0.616f,  0.489f,
  0.997f,  0.513f,  0.064f,
  0.945f,  0.719f,  0.592f,
  0.543f,  0.021f,  0.978f,
  0.279f,  0.317f,  0.505f,
  0.167f,  0.620f,  0.077f,
  0.347f,  0.857f,  0.137f,
  0.055f,  0.953f,  0.042f,
  0.714f,  0.505f,  0.345f,
  0.783f,  0.290f,  0.734f,
  0.722f,  0.645f,  0.174f,
  0.302f,  0.455f,  0.848f,
  0.225f,  0.587f,  0.040f,
  0.517f,  0.713f,  0.338f,
  0.053f,  0.959f,  0.120f,
  0.393f,  0.621f,  0.362f,
  0.673f,  0.211f,  0.457f,
  0.820f,  0.883f,  0.371f,
  0.982f,  0.099f,  0.879f
};

static void parseInput(
    char* filename,
    std::vector<Makena::Vec3>& points
) {
    points.clear();

    std::string   sfilename(filename);
    std::ifstream gs;
    gs.open (sfilename);;

    if (gs.is_open()) {
        std::string   oneLine;
        while (getline(gs, oneLine)) {
            if (oneLine.length()==0)
                continue;
            if (oneLine[0]=='#')
                continue;
            std::stringstream lineStream(oneLine);
            std::string field;
            std::vector<float> p;
            while (std::getline(lineStream, field, ',')) {
                p.push_back(std::atof(field.c_str()));
            }
            Makena::Vec3 v(p[0],p[1],p[2]);
            points.push_back(v);
        }
    }
    gs.close();
}

class GLfloatArrayDelete {public:void operator()(void* x){free(x);}};

using GLfloatArray = std::unique_ptr<GLfloat,GLfloatArrayDelete>;

GLfloatArray makeGLverticesArray(std::vector<Makena::Vec3>& points)
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


GLfloatArray makeGLverticesArray(std::vector<Makena::Vec3>& points, float w)
{
    GLfloatArray ar((GLfloat*)malloc(points.size()*sizeof(GLfloat)*4));
    size_t index =0;
    GLfloat* ap = (GLfloat*)ar.get();
    for (auto& p : points) {
        ap[index++] = p.x();
        ap[index++] = p.y();
        ap[index++] = p.z();
        ap[index++] = w;
    }
    return ar;
}


GLfloatArray makeGLcolorArray(std::vector<Makena::Vec3>& points)
{
    GLfloatArray ar((GLfloat*)malloc(points.size()*sizeof(GLfloat)*3));
    size_t index =0;
    GLfloat* ap = (GLfloat*)ar.get();
    for (int i = 0; i < points.size(); i++) {
        ap[index++] = 0.3;
        ap[index++] = 0.3;
        ap[index++] = 0.3;
    }
    return ar;
}


GLfloatArray makeGLverticesArrayPCAAxes(
    Makena::Mat3x3& pca,
    Makena::Vec3&   spread,
    Makena::Vec3&   mean
) {

    GLfloatArray ar((GLfloat*)malloc(6*sizeof(GLfloat)*3));
    GLfloat* ap = (GLfloat*)ar.get();
    Makena::Vec3 ev1(pca.cell(1,1), pca.cell(2,1), pca.cell(3,1));
    Makena::Vec3 ev2(pca.cell(1,2), pca.cell(2,2), pca.cell(3,2)); 
    Makena::Vec3 ev3(pca.cell(1,3), pca.cell(2,3), pca.cell(3,3)); 
    ev1.normalize();
    ev2.normalize();
    ev3.normalize();
    ev1.scale(sqrt(spread.x()));
    ev2.scale(sqrt(spread.y()));
    ev3.scale(sqrt(spread.z()));

    ev1 += mean;
    ev2 += mean;
    ev3 += mean;
    
    ap[ 0] = ev1.x();
    ap[ 1] = ev1.y();
    ap[ 2] = ev1.z();
    ap[ 3] = mean.x();
    ap[ 4] = mean.y();
    ap[ 5] = mean.z();
    ap[ 6] = ev2.x();
    ap[ 7] = ev2.y();
    ap[ 8] = ev2.z();
    ap[ 9] = mean.x();
    ap[10] = mean.y();
    ap[11] = mean.z();
    ap[12] = ev3.x();
    ap[13] = ev3.y();
    ap[14] = ev3.z();
    ap[15] = mean.x();
    ap[16] = mean.y();
    ap[17] = mean.z();

    return ar;
}


int main( int argc, char* argv[])
{

    std::vector<Makena::Vec3> points;
    parseInput(argv[1],points);
    auto pointsVertices = makeGLverticesArray(points);
    auto pointsColors   = makeGLcolorArray(points);

    Makena::Vec3   spread;
    Makena::Vec3   mean;
    Makena::Mat3x3 pca
            = Makena::findPrincipalComponents(points, spread, mean);

    std::cerr << "Mean: " 
                     << mean.x() << "," << mean.y() << "," << mean.z() << "\n";

    std::cerr << "Spread: " 
               << spread.x() << "," << spread.y() << "," << spread.z() << "\n";

    auto pcaAxesVertices = makeGLverticesArrayPCAAxes(pca, spread, mean);

    TwBar *twBar01;

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

    window = glfwCreateWindow( 1024, 768, "Test Visualiser PCA", NULL, NULL);

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
    TwSetPixelRatio(float(widthInPixel)/(float)width, float(heightInPixel)/(float)height);
    twBar01 = TwNewBar("Control");

    TwDefine(" GLOBAL help='Control Panel' ");

    double distance = 7.0; 
    TwAddVarRW(twBar01, "distance", TW_TYPE_DOUBLE, &distance, 
               " label='Camera Distance' min=2.0 max=20.0 step=0.1 keyIncr=s keyDecr=S help='Rotation speed (turns/second)' ");

    glm::quat orientation;
    TwAddVarRW(twBar01, "Quaternion", TW_TYPE_QUAT4F, &orientation, " showval=true opened=true ");

    TwSetParam(twBar01, NULL, "refresh",  TW_PARAM_CSTRING, 1, "0.1");
    TwSetParam(twBar01, NULL, "position", TW_PARAM_CSTRING, 1, " 10 10 ");

    glfwSetMouseButtonCallback (window, (GLFWmousebuttonfun) TwEventMouseButtonGLFW3);
    glfwSetCursorPosCallback   (window, (GLFWcursorposfun)   TwEventMousePosGLFW3);
    glfwSetScrollCallback      (window, (GLFWscrollfun)      TwEventMouseWheelGLFW3);  
    glfwSetKeyCallback         (window, (GLFWkeyfun)         TwEventKeyGLFW3);
    glfwSetCharCallback        (window, (GLFWcharfun)        TwEventCharGLFW3);


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

    GLuint vertexBuffer;
    glGenBuffers(1, &vertexBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(g_vertex_buffer_data),
        g_vertex_buffer_data,
        GL_STATIC_DRAW
    );

    GLuint colorBuffer;
    glGenBuffers(1, &colorBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, colorBuffer);
    glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(g_color_buffer_data),
        g_color_buffer_data,
        GL_STATIC_DRAW
    );

    GLuint normalBuffer;
    glGenBuffers(1, &normalBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
    glBufferData(
        GL_ARRAY_BUFFER,
        sizeof(g_normal_buffer_data),
        g_normal_buffer_data,
        GL_STATIC_DRAW
    );

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

    GLuint vertexBufferPoints;
    glGenBuffers(1, &vertexBufferPoints);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBufferPoints);
    glBufferData(
        GL_ARRAY_BUFFER,
        points.size()*sizeof(GLfloat)*3,
        pointsVertices.get(),
        GL_STATIC_DRAW
    );

    GLuint colorBufferPoints;
    glGenBuffers(1, &colorBufferPoints);
    glBindBuffer(GL_ARRAY_BUFFER, colorBufferPoints);
    glBufferData(
        GL_ARRAY_BUFFER,
        points.size()*sizeof(GLfloat)*3,
        pointsColors.get(),
        GL_STATIC_DRAW
    );

    GLuint vertexBufferPCAAxes;
    glGenBuffers(1, &vertexBufferPCAAxes);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBufferPCAAxes);
    glBufferData(
        GL_ARRAY_BUFFER,
        6*sizeof(GLfloat)*3,
        pcaAxesVertices.get(),
        GL_STATIC_DRAW
    );

    GLuint SID_PT = loadShaders(
        "src_interactive_tests/shaders/vert_pt.glsl",
        "src_interactive_tests/shaders/frag_pt.glsl" 
    );

    GLuint SID_POL = loadShaders(
        "src_interactive_tests/shaders/vert_pol.glsl",
        "src_interactive_tests/shaders/frag_pol.glsl" 
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

    do {

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


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


        glDisable(GL_BLEND);

        {
            glUseProgram(SID_PT);

            glPointSize(2.0);

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
            glBindBuffer(GL_ARRAY_BUFFER, vertexBufferPoints);
            glVertexAttribPointer(
                         SID_PT_vertexLCS, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, colorBufferPoints);
            glVertexAttribPointer(
                         SID_PT_vertexRGB, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glDrawArrays(GL_POINTS, 0, points.size()*3);

            glDisableVertexAttribArray(1);
            glDisableVertexAttribArray(0);

            glEnableVertexAttribArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, vertexBufferPCAAxes);
            glVertexAttribPointer(
                         SID_PT_vertexLCS, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, colorBufferGCSAxes);
            glVertexAttribPointer(
                         SID_PT_vertexRGB, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glDrawArrays(GL_LINES, 0, 6*3);

            glDisableVertexAttribArray(1);
            glDisableVertexAttribArray(0);

        }

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


        {
            glUseProgram(SID_POL);

            glUniformMatrix4fv(SID_POL_MVP, 1, GL_FALSE, &MVP[0][0] );
            glUniformMatrix4fv(SID_POL_M,   1, GL_FALSE, &M  [0][0] );
            glUniformMatrix4fv(SID_POL_V,   1, GL_FALSE, &V  [0][0] );
            glUniform3f(SID_POL_lightWCS, lightPos.x, lightPos.y, lightPos.z);

            glEnableVertexAttribArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
            glVertexAttribPointer(
                         SID_POL_vertexLCS, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, colorBuffer);
            glVertexAttribPointer(
                        SID_POL_vertexRGBA, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

            glEnableVertexAttribArray(2);
            glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
            glVertexAttribPointer(
                         SID_POL_normalLCS, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

//            glDrawArrays(GL_TRIANGLES, 0, 12*3);

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

    glDeleteBuffers(1, &vertexBuffer);
    glDeleteBuffers(1, &colorBuffer);
    glDeleteBuffers(1, &normalBuffer);
    glDeleteVertexArrays(1, &VertexArrayID);

    glDeleteProgram(SID_POL);
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


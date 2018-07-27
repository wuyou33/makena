#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
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

#include "constraint_manager.hpp"
#include "contact_manager.hpp"

class GlfwManager;
class TwManager;

/** @file src_interactive_tests/test_visualizer_contacts.cpp
 *
 *  @brief visual test rig for Makena::ContactManager
 *
 *  Usage: test_visualizer_contacts
 *
 *  @dependencies
 *    OpenGL 3.3 or later
 *    GLFW3
 *    GLEW
 *    GLM
 *    AntTweakbar 1.16 custom (updated for GLFW3 and Macbook Retina)
 */


/** @brief subcomponent class definitions
 *         - TwManager
 *         - GlfwManager
 *         - SharderManager
 *         - GLBufferManager
 *
 */


class TwCallable {
  public:
    virtual void TwUpdate(TwManager& tw) = 0;
};


class TwManager {

  public:

    TwManager(GlfwManager& glfw);

    virtual ~TwManager();

    virtual void init();

    void initCommon();

    void setCallback(TwCallable* c);

    void draw();

    void terminate();

    void updateVers();

    static void TwEventMouseButtonGLFW3(
                        GLFWwindow* window, int button, int action, int mods);

    static void TwEventMousePosGLFW3(
                        GLFWwindow* window, double xpos, double ypos);

    static void TwEventMouseWheelGLFW3(
                        GLFWwindow* window, double xoffset, double yoffset);

    static void TwEventKeyGLFW3(
              GLFWwindow* window, int key, int scancode, int action, int mods);

    static void TwEventCharGLFW3(GLFWwindow* window, int codepoint);

    GlfwManager& mGlfw;

    TwBar*       mTwMainBar;

    double       mViewDistance;
    double       mViewVerticalPosition;
    double       mViewHorizontalPosition;
    glm::quat    mViewOrientation;

    TwCallable*  mCallable;
};


class GlfwCallable {
  public:
    virtual void GlfwStep() = 0;
    virtual void GlfwKeyCallback(char ch) = 0;
};


class GlfwManager  {

  public:

    GlfwManager(int requestedWindowWidth,int requestedWindowHeight);

    bool initGLFW();

    void configGLFW();

    void setCallback(GlfwCallable* c);

    int windowWidth()         const { return mWindowWidth; }
    int windowHeight()        const { return mWindowHeight; }
    int windowWidthInPixel()  const { return mWindowWidthInPixel; }
    int windowHeightInPixel() const { return mWindowHeightInPixel; }

    void update();

    bool shouldExit();

    void terminate();

  private:

    static void callbackWindowSize(GLFWwindow *window, int w, int h);
    static void callbackCursorPos (GLFWwindow *window, double x, double y);

    GLFWwindow*   mWindow;

    int           mRequestedWindowWidth;
    int           mRequestedWindowHeight;
    int           mWindowWidth;
    int           mWindowHeight;
    int           mWindowWidthInPixel;
    int           mWindowHeightInPixel;

    bool          mSpacePressed;
    bool          mAPressed;
    bool          mBPressed;
    bool          mCPressed;
    bool          mDPressed;
    bool          mEPressed;
    bool          mFPressed;
    bool          mGPressed;
    bool          mHPressed;
    bool          mIPressed;
    bool          mJPressed;
    bool          mKPressed;
    bool          mLPressed;
    bool          mMPressed;
    bool          mNPressed;
    bool          mOPressed;
    bool          mPPressed;
    bool          mQPressed;
    bool          mRPressed;
    bool          mSPressed;
    bool          mTPressed;
    bool          mUPressed;
    bool          mVPressed;
    bool          mWPressed;
    bool          mXPressed;
    bool          mYPressed;
    bool          mZPressed;
    GlfwCallable* mCallable;

};


class ShaderManager {

  public:
    ShaderManager (TwManager& TW);
    ~ShaderManager ();

    void loadShaders();

    void setViewGeometry();

    void loadParamsToGL(const GLuint shaderType);

    void unloadShaders();



    GLuint SID_PT;
    GLuint SID_POL;

    GLuint SID_POL_MVP;
    GLuint SID_POL_M;
    GLuint SID_POL_V;
    GLuint SID_POL_lightWCS;

    GLuint SID_POL_vertexLCS;
    GLuint SID_POL_vertexRGBA;
    GLuint SID_POL_normalLCS;

    GLuint SID_PT_MVP;
    GLuint SID_PT_M;
    GLuint SID_PT_lightWCS;

    GLuint SID_PT_vertexLCS;
    GLuint SID_PT_vertexRGB;

    glm::mat4 M;
    glm::mat4 V;
    glm::mat4 MVP;
    glm::vec3 lightPos;

  private:
    TwManager& mTW;

    GLuint loadShadersForOneType(const char * v_path,const char * f_path);

};


class GLBufferManager {

  public:
    GLBufferManager(ShaderManager& mShaders);
    

    class GLfloatArrayDelete {public:void operator()(void* x){free(x);}};

    using GLfloatArray = std::unique_ptr<GLfloat,GLfloatArrayDelete>;

    void genBuffers(long capacity);

    void releaseBuffers();

    void resetOneStep();

    void addTriangles(
        const std::vector<Makena::Vec3>& p,
        const std::vector<Makena::Vec3>& c,
        const std::vector<float>&        a,
        const std::vector<Makena::Vec3>& n
    );

    void addLines(
        const std::vector<Makena::Vec3>& p,
        const std::vector<Makena::Vec3>& c
    );

    void addPoints(
        const std::vector<Makena::Vec3>& p,
        const std::vector<Makena::Vec3>& c
    );

    void draw();

  private:
    void renderTriangles();
    void renderLines();
    void renderPoints();

    GLfloatArray makeGLbufferArray3D(const std::vector<Makena::Vec3>& points);
    GLfloatArray makeGLbufferArrayRGBA(
        const std::vector<Makena::Vec3>& colors,
        const std::vector<float>&        alphas
    );

    ShaderManager&            mShaders;
    GLuint                    mIdVertexArray;
    GLuint                    mIdPointsVertices;
    GLuint                    mIdPointsColors;
    GLuint                    mIdLinesVertices;
    GLuint                    mIdLinesColors;
    GLuint                    mIdTrianglesVertices;
    GLuint                    mIdTrianglesColors;
    GLuint                    mIdTrianglesNormals;

    std::vector<Makena::Vec3> mPointsTriangles;
    std::vector<Makena::Vec3> mColorsTriangles;
    std::vector<float>        mAlphasTriangles;
    std::vector<Makena::Vec3> mNormalsTriangles;
    std::vector<Makena::Vec3> mPointsLines;
    std::vector<Makena::Vec3> mColorsLines;
    std::vector<Makena::Vec3> mPointsPoints;
    std::vector<Makena::Vec3> mColorsPoints;


};


/***************************************************************************/
/**                 @brief TwManager implementation BEGIN                  */
/***************************************************************************/


TwManager::TwManager(GlfwManager& glfw):mGlfw(glfw),mCallable(nullptr){;}


void TwManager::initCommon()
{

    TwInit(TW_OPENGL_CORE, NULL);

    TwWindowSize(mGlfw.windowWidthInPixel(), mGlfw.windowHeightInPixel());

    TwSetPixelRatio(
            (float)mGlfw.windowWidthInPixel()  /(float)mGlfw.windowWidth(),
            (float)mGlfw.windowHeightInPixel() /(float)mGlfw.windowHeight() );

    mTwMainBar = TwNewBar("Main Control");

    TwDefine(" GLOBAL help='Control Panel' ");


    mViewDistance            = 10.0;
    mViewVerticalPosition    = 0.0;
    mViewHorizontalPosition  = 0.0;;

    TwAddVarRW(
        mTwMainBar,
        "view orientation",
        TW_TYPE_QUAT4F,
        &mViewOrientation, 
        " showval=true opened=true ");

    TwAddVarRW(
        mTwMainBar,
        "view distance",
        TW_TYPE_DOUBLE,
        &mViewDistance, 
        " label='View Distance' min=2.0 max=50.0 step=0.1 keyIncr=s "
        "keyDecr=S help='Rotation speed (turns/second)' ");

    TwAddVarRW(
        mTwMainBar,
        "view vertical position",
        TW_TYPE_DOUBLE,
        &mViewVerticalPosition,
        " label='View Vertical Position' min=-100.0 max=100.0 step=0.1 "
        "keyIncr=s keyDecr=S help='Rotation speed (turns/second)' ");

    TwAddVarRW(
        mTwMainBar,
        "view horizontal position",
        TW_TYPE_DOUBLE,
        &mViewHorizontalPosition,
        " label='View Horizontal Position' min=-100.0 max=100.0 step=0.1 "
        "keyIncr=s keyDecr=S help='Rotation speed (turns/second)' ");

}

void TwManager::init()
{
    initCommon();
}


void TwManager::draw(){ TwDraw(); }


void TwManager::terminate(){ TwTerminate(); }


TwManager::~TwManager(){;}


void TwManager::setCallback(TwCallable* c)
{
    mCallable = c;
}


void TwManager::updateVers()
{
    if (mCallable != nullptr) {
        mCallable->TwUpdate(*this);
    }
}



void TwManager::TwEventMouseButtonGLFW3(
    GLFWwindow* window,
    int         button,
    int         action,
    int         mods
) {
    TwEventMouseButtonGLFW(button, action);
}


void TwManager::TwEventMousePosGLFW3(
    GLFWwindow* window,
    double      xpos,
    double      ypos
) {
    TwMouseMotion(int(xpos), int(ypos));
}


void TwManager::TwEventMouseWheelGLFW3(
    GLFWwindow* window,
    double      xoffset,
    double      yoffset
) {
    TwEventMouseWheelGLFW(yoffset);
}


void TwManager::TwEventKeyGLFW3(
    GLFWwindow* window,
    int         key,
    int         scancode,
    int         action,
    int         mods
) {
    TwEventKeyGLFW(key, action);
}


void TwManager::TwEventCharGLFW3(GLFWwindow* window, int codepoint)
{
    TwEventCharGLFW(codepoint, GLFW_PRESS);
}


/***************************************************************************/
/**                 @brief TwManager implementation END                    */
/***************************************************************************/


/***************************************************************************/
/**                 @brief GlfwManager implementation BEGIN                */
/***************************************************************************/


GlfwManager::GlfwManager(
    int requestedWindowWidth,
    int requestedWindowHeight
):
    mRequestedWindowWidth(requestedWindowWidth),
    mRequestedWindowHeight(requestedWindowHeight),
    mSpacePressed(false),
    mAPressed(false),
    mBPressed(false),
    mCPressed(false),
    mDPressed(false),
    mEPressed(false),
    mFPressed(false),
    mGPressed(false),
    mHPressed(false),
    mIPressed(false),
    mJPressed(false),
    mKPressed(false),
    mLPressed(false),
    mMPressed(false),
    mNPressed(false),
    mOPressed(false),
    mPPressed(false),
    mQPressed(false),
    mRPressed(false),
    mSPressed(false),
    mTPressed(false),
    mUPressed(false),
    mVPressed(false),
    mWPressed(false),
    mXPressed(false),
    mYPressed(false),
    mZPressed(false),
    mCallable(nullptr){;}


bool GlfwManager::initGLFW()
{

    // Initialise GLFW
    if (!glfwInit()) {
        std::cerr << "glfwInit() failed.\n";
        return false;
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    mWindow = glfwCreateWindow(
                 mRequestedWindowWidth, mRequestedWindowHeight,
                 "Test Visualiser Convex Hull", NULL, NULL);

    if ( mWindow == NULL ) {
        std::cerr << "glfwCreateWindow failed\n";
        glfwTerminate();
        return false;
    }

    glfwMakeContextCurrent(mWindow);
    
    glfwGetFramebufferSize(
                      mWindow, &mWindowWidthInPixel, &mWindowHeightInPixel);

    glfwSetWindowSizeCallback(mWindow, callbackWindowSize);

    glfwGetWindowSize(mWindow, &mWindowWidth, &mWindowHeight);

    glewExperimental = true;
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW\n";;
        terminate();
        return false;
    }

    return true;
}

void GlfwManager::configGLFW()
{

    glfwSetMouseButtonCallback (mWindow, (GLFWmousebuttonfun) 
                                        TwManager::TwEventMouseButtonGLFW3);
    glfwSetCursorPosCallback   (mWindow, (GLFWcursorposfun)
                                        TwManager::TwEventMousePosGLFW3);
    glfwSetScrollCallback      (mWindow, (GLFWscrollfun) 
                                        TwManager::TwEventMouseWheelGLFW3);
    glfwSetKeyCallback         (mWindow, (GLFWkeyfun) 
                                        TwManager::TwEventKeyGLFW3);
    glfwSetCharCallback        (mWindow, (GLFWcharfun)
                                        TwManager::TwEventCharGLFW3);

    glfwSetInputMode(mWindow, GLFW_STICKY_KEYS, GL_TRUE);
    glfwSetInputMode(mWindow, GLFW_CURSOR,      GLFW_CURSOR_NORMAL);

    glfwPollEvents();

    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glDepthFunc(GL_LESS); 

}


void GlfwManager::setCallback(GlfwCallable* c)
{
    mCallable = c;
}

void GlfwManager::update()
{
    glfwSwapBuffers(mWindow);
    glfwPollEvents();

    if (glfwGetKey(mWindow,GLFW_KEY_SPACE ) == GLFW_PRESS && !mSpacePressed) {
        mSpacePressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwStep();
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_SPACE)!=GLFW_PRESS &&mSpacePressed) {
        mSpacePressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_A)==GLFW_PRESS &&!mAPressed) {
        mAPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('a');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_A)!=GLFW_PRESS &&mAPressed) {
        mAPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_B)==GLFW_PRESS &&!mBPressed) {
        mBPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('b');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_B)!=GLFW_PRESS &&mBPressed) {
        mBPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_C)==GLFW_PRESS &&!mCPressed) {
        mCPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('c');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_C)!=GLFW_PRESS &&mCPressed) {
        mCPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_D)==GLFW_PRESS &&!mDPressed) {
        mDPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('d');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_D)!=GLFW_PRESS &&mDPressed) {
        mDPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_E)==GLFW_PRESS &&!mEPressed) {
        mEPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('d');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_E)!=GLFW_PRESS &&mEPressed) {
        mEPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_F)==GLFW_PRESS &&!mFPressed) {
        mFPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('f');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_F)!=GLFW_PRESS &&mFPressed) {
        mFPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_G)==GLFW_PRESS &&!mGPressed) {
        mGPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('g');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_G)!=GLFW_PRESS &&mGPressed) {
        mGPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_H)==GLFW_PRESS &&!mHPressed) {
        mHPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('h');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_H)!=GLFW_PRESS &&mHPressed) {
        mHPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_I)==GLFW_PRESS &&!mIPressed) {
        mIPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('i');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_I)!=GLFW_PRESS &&mIPressed) {
        mIPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_J)==GLFW_PRESS &&!mJPressed) {
        mJPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('j');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_J)!=GLFW_PRESS &&mJPressed) {
        mJPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_K)==GLFW_PRESS &&!mKPressed) {
        mKPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('k');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_K)!=GLFW_PRESS &&mKPressed) {
        mKPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_L)==GLFW_PRESS &&!mLPressed) {
        mLPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('l');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_L)!=GLFW_PRESS &&mLPressed) {
        mLPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_M)==GLFW_PRESS &&!mMPressed) {
        mMPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('m');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_M)!=GLFW_PRESS &&mMPressed) {
        mMPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_N)==GLFW_PRESS &&!mNPressed) {
        mNPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('n');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_N)!=GLFW_PRESS &&mNPressed) {
        mNPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_O)==GLFW_PRESS &&!mOPressed) {
        mOPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('o');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_O)!=GLFW_PRESS &&mOPressed) {
        mOPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_P)==GLFW_PRESS &&!mPPressed) {
        mPPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('p');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_P)!=GLFW_PRESS &&mPPressed) {
        mPPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_Q)==GLFW_PRESS &&!mQPressed) {
        mQPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('q');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_Q)!=GLFW_PRESS &&mQPressed) {
        mQPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_R)==GLFW_PRESS &&!mRPressed) {
        mRPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('r');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_R)!=GLFW_PRESS &&mRPressed) {
        mRPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_S)==GLFW_PRESS &&!mSPressed) {
        mSPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('s');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_S)!=GLFW_PRESS &&mSPressed) {
        mSPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_T)==GLFW_PRESS &&!mTPressed) {
        mTPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('t');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_T)!=GLFW_PRESS &&mTPressed) {
        mTPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_U)==GLFW_PRESS &&!mUPressed) {
        mUPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('u');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_U)!=GLFW_PRESS &&mUPressed) {
        mUPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_V)==GLFW_PRESS &&!mVPressed) {
        mVPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('v');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_V)!=GLFW_PRESS &&mVPressed) {
        mVPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_W)==GLFW_PRESS &&!mWPressed) {
        mWPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('w');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_W)!=GLFW_PRESS &&mWPressed) {
        mWPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_X)==GLFW_PRESS &&!mXPressed) {
        mXPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('x');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_X)!=GLFW_PRESS &&mXPressed) {
        mXPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_Y)==GLFW_PRESS &&!mYPressed) {
        mYPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('y');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_Y)!=GLFW_PRESS &&mYPressed) {
        mYPressed = false;
    }

    if (glfwGetKey(mWindow,GLFW_KEY_Z)==GLFW_PRESS &&!mZPressed) {
        mZPressed = true;
        if (mCallable!=nullptr) {
            mCallable->GlfwKeyCallback('z');
        }
    }
    else if (glfwGetKey(mWindow,GLFW_KEY_Z)!=GLFW_PRESS &&mZPressed) {
        mZPressed = false;
    }
}


bool GlfwManager::shouldExit()
{
    return glfwGetKey(mWindow, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
           glfwWindowShouldClose(mWindow)        == 0          ;
}


void GlfwManager::terminate()
{
    glfwTerminate();
}


void GlfwManager::callbackWindowSize(GLFWwindow *window, int w, int h)
{

    int widthInPixel, heightInPixel;

    glfwGetFramebufferSize(window, &widthInPixel, &heightInPixel);

    glViewport( 0, 0, (GLsizei)widthInPixel, (GLsizei)heightInPixel);

    TwWindowSize(widthInPixel, heightInPixel);

}


void GlfwManager::callbackCursorPos(GLFWwindow *window, double x, double y)
{
    //std::cerr << "CB cursor pos: " << x << "," << y << "\n";
}


/***************************************************************************/
/**                 @brief GlfwManager implementation END                  */
/***************************************************************************/


/***************************************************************************/
/**                 @brief ShaderManager implementation BEGIN              */
/***************************************************************************/

ShaderManager::ShaderManager(TwManager& TW):mTW(TW){;}

ShaderManager::~ShaderManager () {;}

void ShaderManager::loadShaders()
{

    SID_PT = loadShadersForOneType(
        "src_interactive_tests/shaders/vert_pt.glsl",
        "src_interactive_tests/shaders/frag_pt.glsl" 
    );

    SID_POL = loadShadersForOneType(
        "src_interactive_tests/shaders/vert_pol.glsl",
        "src_interactive_tests/shaders/frag_pol.glsl" 
    );

    SID_POL_MVP        = glGetUniformLocation(SID_POL, "MVP");
    SID_POL_M          = glGetUniformLocation(SID_POL, "M");
    SID_POL_V          = glGetUniformLocation(SID_POL, "V");
    SID_POL_lightWCS   = glGetUniformLocation(SID_POL, "lightWCS");

    SID_POL_vertexLCS  = glGetAttribLocation(SID_POL, "vertexLCS");
    SID_POL_vertexRGBA = glGetAttribLocation(SID_POL, "vertexRGBA");
    SID_POL_normalLCS  = glGetAttribLocation(SID_POL, "normalLCS");

    SID_PT_MVP         = glGetUniformLocation(SID_PT, "MVP");
    SID_PT_M           = glGetUniformLocation(SID_PT, "M");
    SID_PT_lightWCS    = glGetUniformLocation(SID_PT, "lightWCS");

    SID_PT_vertexLCS   = glGetAttribLocation(SID_PT, "vertexLCS");
    SID_PT_vertexRGB   = glGetAttribLocation(SID_PT, "vertexRGB");
}


GLuint ShaderManager::loadShadersForOneType(
    const char * v_path,const char * f_path
) {
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

void ShaderManager::setViewGeometry()
{
    glm::vec3 scaleFactor(1.0f, 1.0f, 1.0f);
    glm::mat4 Mscale = glm::scale( glm::mat4(), scaleFactor );
                                   
    glm::quat Qrot(0.0f, 1.0f, 0.0f, 0.0f);
    glm::mat4 Mrot   = glm::mat4_cast(mTW.mViewOrientation);

    glm::mat4 Mtrans = glm::translate(glm::mat4(), 
                                      glm::vec3(
                                          mTW.mViewHorizontalPosition,
                                          mTW.mViewVerticalPosition,
                                          0.0f
                                      )
                                     );
    M = Mtrans * Mrot * Mscale;

    V = glm::lookAt(glm::vec3(0.0f, 0.0f, mTW.mViewDistance ),// Cam pos
                    glm::vec3( 0.0f, 0.0f, 0.0f ),            // and looks here
                    glm::vec3( 0.0f, 1.0f, 0.0f ) );          // Head is up

    glm::mat4 P      = glm::perspective(
                               glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.0f);
 
    MVP              = P * V * M;

    lightPos         = glm::vec3(4.0, 4.0, 4.0);
}


void ShaderManager::loadParamsToGL(const GLuint shaderType)
{
    if (shaderType == SID_PT) {

        glUniformMatrix4fv(SID_PT_MVP, 1, GL_FALSE, &MVP[0][0] );
        glUniformMatrix4fv(SID_PT_M,   1, GL_FALSE, &M  [0][0] );
        glUniform3f(SID_PT_lightWCS, lightPos.x, lightPos.y, lightPos.z);

    }                  
    else if (shaderType == SID_POL) {

        glUniformMatrix4fv(SID_POL_MVP,1, GL_FALSE, &MVP[0][0] );
        glUniformMatrix4fv(SID_POL_M,  1, GL_FALSE, &M  [0][0] );
        glUniformMatrix4fv(SID_POL_V,  1, GL_FALSE, &V  [0][0] );
        glUniform3f(SID_POL_lightWCS, lightPos.x, lightPos.y, lightPos.z);

    }
}


void ShaderManager::unloadShaders() {    

    glDeleteProgram(SID_POL);
    glDeleteProgram(SID_PT);

}


/***************************************************************************/
/**                 @brief ShaderManager implementation END                */
/***************************************************************************/


/***************************************************************************/
/**               @brief GLBufferManagerimplementation BEGIN               */
/***************************************************************************/


GLBufferManager::GLBufferManager(ShaderManager& shaders):mShaders(shaders)
{;}


void GLBufferManager::genBuffers(long capacity) {

    glGenVertexArrays(1, &mIdVertexArray);

    glGenBuffers(1, &mIdPointsVertices);
    glGenBuffers(1, &mIdPointsColors);
    glGenBuffers(1, &mIdLinesVertices);
    glGenBuffers(1, &mIdLinesColors);
    glGenBuffers(1, &mIdTrianglesVertices);
    glGenBuffers(1, &mIdTrianglesColors);
    glGenBuffers(1, &mIdTrianglesNormals);

}


void GLBufferManager::releaseBuffers() {

    glDeleteBuffers(1, &mIdTrianglesVertices);
    glDeleteBuffers(1, &mIdTrianglesColors);
    glDeleteBuffers(1, &mIdTrianglesNormals);
    glDeleteBuffers(1, &mIdLinesVertices);
    glDeleteBuffers(1, &mIdLinesColors);
    glDeleteBuffers(1, &mIdPointsVertices);
    glDeleteBuffers(1, &mIdPointsColors);

    glDeleteVertexArrays(1, &mIdVertexArray);

}


void GLBufferManager::resetOneStep()
{
    mPointsTriangles.clear();
    mColorsTriangles.clear();
    mAlphasTriangles.clear();
    mNormalsTriangles.clear();
    mPointsLines.clear();
    mColorsLines.clear();
    mPointsPoints.clear();
    mColorsPoints.clear();
}


void GLBufferManager::addTriangles(
    const std::vector<Makena::Vec3>& p,
    const std::vector<Makena::Vec3>& c,
    const std::vector<float>&        a,
    const std::vector<Makena::Vec3>& n
) {

    mPointsTriangles.insert(mPointsTriangles.begin(),   p.begin(), p.end());
    mColorsTriangles.insert(mColorsTriangles.begin(),   c.begin(), c.end());
    mAlphasTriangles.insert(mAlphasTriangles.begin(),   a.begin(), a.end());
    mNormalsTriangles.insert(mNormalsTriangles.begin(), n.begin(), n.end());

}


void GLBufferManager::addLines(
    const std::vector<Makena::Vec3>& p,
    const std::vector<Makena::Vec3>& c
) {

    mPointsLines.insert(mPointsLines.begin(),   p.begin(), p.end());
    mColorsLines.insert(mColorsLines.begin(),   c.begin(), c.end());
}


void GLBufferManager::addPoints(
    const std::vector<Makena::Vec3>& p,
    const std::vector<Makena::Vec3>& c
) {

    mPointsPoints.insert(mPointsPoints.begin(),   p.begin(), p.end());
    mColorsPoints.insert(mColorsPoints.begin(),   c.begin(), c.end());
}


void GLBufferManager::draw()
{

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (!mPointsTriangles.empty()) {
        renderTriangles()   ;
    }

    if (!mPointsLines.empty()) {
        renderLines()   ;
    }

    if (!mPointsPoints.empty()) {
        renderPoints()   ;
    }

}


GLBufferManager::GLfloatArray GLBufferManager::makeGLbufferArray3D(
    const std::vector<Makena::Vec3>& points
) {
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


GLBufferManager::GLfloatArray GLBufferManager::makeGLbufferArrayRGBA(
    const std::vector<Makena::Vec3>& colors,
    const std::vector<float>&        alphas
) {
    GLfloatArray ar((GLfloat*)malloc(colors.size()*sizeof(GLfloat)*4));
    size_t index =0;
    size_t alphaIndex = 0;
    GLfloat* ap = (GLfloat*)ar.get();
    for (auto& p : colors) {
        ap[index++] = p.x();
        ap[index++] = p.y();
        ap[index++] = p.z();
        ap[index++] = alphas[alphaIndex++];
    }
    return ar;
}


void GLBufferManager::renderTriangles()
{

    auto arrTrianglesVertices = makeGLbufferArray3D  ( mPointsTriangles  );
    auto arrTrianglesColors   = makeGLbufferArrayRGBA( mColorsTriangles, 
                                                       mAlphasTriangles  );
    auto arrTrianglesNormals  = makeGLbufferArray3D  ( mNormalsTriangles );

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glUseProgram(mShaders.SID_POL);

    mShaders.loadParamsToGL(mShaders.SID_POL);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, mIdTrianglesVertices);
    glBufferData(
        GL_ARRAY_BUFFER,
        mPointsTriangles.size()*sizeof(GLfloat)*3,
        arrTrianglesVertices.get(),
        GL_DYNAMIC_DRAW
    );
    glVertexAttribPointer(
            mShaders.SID_POL_vertexLCS, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, mIdTrianglesColors);
    glBufferData(
        GL_ARRAY_BUFFER,
        mColorsTriangles.size()*sizeof(GLfloat)*4,
        arrTrianglesColors.get(),
        GL_DYNAMIC_DRAW
    );

    glVertexAttribPointer(
               mShaders.SID_POL_vertexRGBA, 4, GL_FLOAT, GL_FALSE, 0,(void*)0);

    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, mIdTrianglesNormals);
    glBufferData(
        GL_ARRAY_BUFFER,
        mNormalsTriangles.size()*sizeof(GLfloat)*3,
        arrTrianglesNormals.get(),
        GL_DYNAMIC_DRAW
    );

    glVertexAttribPointer(
                mShaders.SID_POL_normalLCS, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);

    glDrawArrays(GL_TRIANGLES, 0, mPointsTriangles.size());

    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(0);

}


void GLBufferManager::renderLines()
{

    auto arrLinesVertices = makeGLbufferArray3D( mPointsLines );
    auto arrLinesColors   = makeGLbufferArray3D( mColorsLines );

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    glDisable(GL_BLEND);
    glUseProgram(mShaders.SID_PT);
    glLineWidth(10.0);

    mShaders.loadParamsToGL(mShaders.SID_PT);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, mIdPointsVertices);
    glBufferData(
        GL_ARRAY_BUFFER,
        mPointsLines.size()*sizeof(GLfloat)*3,
        arrLinesVertices.get(),
        GL_DYNAMIC_DRAW
    );

    glVertexAttribPointer(
                 mShaders.SID_PT_vertexLCS, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);
    
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, mIdPointsColors);
    glBufferData(
        GL_ARRAY_BUFFER,
        mColorsLines.size()*sizeof(GLfloat)*3,
        arrLinesColors.get(),
        GL_DYNAMIC_DRAW
    );

    glVertexAttribPointer(
                 mShaders.SID_PT_vertexRGB, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);
                         
    glDrawArrays(GL_LINES, 0, mPointsLines.size());

    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(0);
}


void GLBufferManager::renderPoints()
{
    auto arrPointsVertices = makeGLbufferArray3D( mPointsPoints );
    auto arrPointsColors   = makeGLbufferArray3D( mColorsPoints );

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glDisable(GL_BLEND);
    glUseProgram(mShaders.SID_PT);
    glPointSize(10.0);

    mShaders.loadParamsToGL(mShaders.SID_PT);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, mIdPointsVertices);
    glBufferData(
        GL_ARRAY_BUFFER,
        mPointsPoints.size()*sizeof(GLfloat)*3,
        arrPointsVertices.get(),
        GL_DYNAMIC_DRAW
    );

    glVertexAttribPointer(
                 mShaders.SID_PT_vertexLCS, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);
    
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, mIdPointsColors);
    glBufferData(
        GL_ARRAY_BUFFER,
        mColorsPoints.size()*sizeof(GLfloat)*3,
        arrPointsColors.get(),
        GL_DYNAMIC_DRAW
    );

    glVertexAttribPointer(
                 mShaders.SID_PT_vertexRGB, 3, GL_FLOAT, GL_FALSE, 0,(void*)0);
                         
    glDrawArrays(GL_POINTS, 0, mPointsPoints.size());

    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(0);
}


/***************************************************************************/
/**               @brief GLBufferManagerimplementation END                 */
/***************************************************************************/


/***************************************************************************/
/**               @brief GCSAxesLines BEGIN                                */
/***************************************************************************/


class GCSAxesLines {

  public:
    GCSAxesLines(float markerLength, float totalLength);
    ~GCSAxesLines();

    std::vector<Makena::Vec3> mVertices;
    std::vector<Makena::Vec3> mColors;

};



GCSAxesLines::GCSAxesLines(float markerLength, float totalLength)
{
    Makena::Vec3 ori(0.0, 0.0, 0.0);

    Makena::Vec3 x1(    markerLength, 0.0, 0.0);
    Makena::Vec3 x2( 0.5*totalLength, 0.0, 0.0);
    Makena::Vec3 x3(-0.0*totalLength, 0.0, 0.0);

    Makena::Vec3 y1(0.0,     markerLength, 0.0);
    Makena::Vec3 y2(0.0,  0.5*totalLength, 0.0);
    Makena::Vec3 y3(0.0, -0.0*totalLength, 0.0);

    Makena::Vec3 z1(0.0, 0.0,     markerLength);
    Makena::Vec3 z2(0.0, 0.0,  0.5*totalLength);
    Makena::Vec3 z3(0.0, 0.0, -0.0*totalLength);

    Makena::Vec3 red_high  (1.0, 0.0, 0.0);
    Makena::Vec3 red_low   (0.2, 0.0, 0.0);
    Makena::Vec3 green_high(0.0, 1.0, 0.0);
    Makena::Vec3 green_low (0.0, 0.2, 0.0);
    Makena::Vec3 blue_high (0.0, 0.0, 1.0);
    Makena::Vec3 blue_low  (0.0, 0.0, 0.2);

    mVertices.push_back(ori);
    mVertices.push_back(x1);
    mColors.push_back(red_high);
    mColors.push_back(red_high);

    mVertices.push_back(x1);
    mVertices.push_back(x2);
    mColors.push_back(red_low);
    mColors.push_back(red_low);

    mVertices.push_back(ori);
    mVertices.push_back(x3);
    mColors.push_back(red_low);
    mColors.push_back(red_low);

    mVertices.push_back(ori);
    mVertices.push_back(y1);
    mColors.push_back(green_high);
    mColors.push_back(green_high);

    mVertices.push_back(y1);
    mVertices.push_back(y2);
    mColors.push_back(green_low);
    mColors.push_back(green_low);

    mVertices.push_back(ori);
    mVertices.push_back(y3);
    mColors.push_back(green_low);
    mColors.push_back(green_low);

    mVertices.push_back(ori);
    mVertices.push_back(z1);
    mColors.push_back(blue_high);
    mColors.push_back(blue_high);

    mVertices.push_back(z1);
    mVertices.push_back(z2);
    mColors.push_back(blue_low);
    mColors.push_back(blue_low);

    mVertices.push_back(ori);
    mVertices.push_back(z3);
    mColors.push_back(blue_low);
    mColors.push_back(blue_low);
}

GCSAxesLines::~GCSAxesLines(){;}


/***************************************************************************/
/**               @brief GCSAxesLines END                                  */
/***************************************************************************/


/***************************************************************************/
/**                 @brief TWManagerAppSpecific BEGIN                      */
/***************************************************************************/


class TwManagerAppSpecific : public TwManager {

  public:

    TwManagerAppSpecific(GlfwManager& glfw);

    virtual ~TwManagerAppSpecific();

    virtual void init();

    glm::quat      mGravityOri;

};

TwManagerAppSpecific::TwManagerAppSpecific(GlfwManager&glfw):
    TwManager(glfw){;}




TwManagerAppSpecific::~TwManagerAppSpecific() {;}

void TwManagerAppSpecific::init()
{

    initCommon();

    TwAddVarRW(mTwMainBar, "Gravity Direction", TW_TYPE_QUAT4F, &mGravityOri,
                                               " showval=true opened=true ");
    TwSetParam(mTwMainBar, NULL, "refresh",  TW_PARAM_CSTRING, 1, "0.1");
    TwSetParam(mTwMainBar, NULL, "position", TW_PARAM_CSTRING, 1, " 10 10 ");
}


/***************************************************************************/
/**                   @brief TWManagerAppSpecific END                      */
/***************************************************************************/


/***************************************************************************/
/**                   @brief TestRig BEGIN                                 */
/***************************************************************************/


class TestRig : public GlfwCallable, public TwCallable {

  public:

    TestRig();
    ~TestRig();

    void init();
    void step(const double deltaT);
    void render(GLBufferManager& buf);
    void term();

    void GlfwStep();

    void GlfwKeyCallback(char ch);

    void GlfwFindClosestP();

    void TwUpdate(TwManager& tw);

  private:

    void constructBase();
    void constructFixed01();
    void constructFixed02();

    Makena::Vec3              mG;
    Makena::ConstraintManager mConsM;
    Makena::ContactManager    mContM;
    Makena::ConvexRigidBody*  mBase;
    Makena::ConvexRigidBody*  mFixed01;
    Makena::ConvexRigidBody*  mFixed02;
    Makena::ConvexRigidBody*  mBody01;
    Makena::ConvexRigidBody*  mBody02;
    Makena::ConvexRigidBody*  mBody03;
    Makena::ConvexRigidBody*  mBody04;
    Makena::ConvexRigidBody*  mBody05;
    bool                      mGo;

};


TestRig::TestRig():
    mContM(0.90,   // KCorr,
        0.01,    // velocityThreshold,
        6,       // maxNumPointsPerContact,
        1.0e-06, // epsilonZero,
        0.1,     // epsilonAngle,
        1.0e-06, // epsilonZeroGJK,
        1.0,     // scalingGJK,
        100,     // maxNumIter,
        8,       // maxNumCycles,
        std::cerr// logStream
    ),
    mBase(nullptr),
    mBody01(nullptr),
    mBody02(nullptr),
    mBody03(nullptr),
    mBody04(nullptr),
    mBody05(nullptr),
    mGo(false)
    {;}

TestRig::~TestRig(){;}

void TestRig::GlfwStep(){;}
void TestRig::GlfwKeyCallback(char ch){ mGo=true; }


void TestRig::TwUpdate(TwManager& tw)
{
    auto& tw2 = dynamic_cast<TwManagerAppSpecific&>(tw);

    Makena::Quaternion q1( tw2.mGravityOri.w,
                           tw2.mGravityOri.x,
                           tw2.mGravityOri.y,
                           tw2.mGravityOri.z  );
    q1.normalize();
    Makena::Vec3 dir(0.0, -1.0, 0.0);
//    Makena::Vec3 dir(1.0, -1.0, -1.0);
    dir.normalize();
    mG = q1.rotate(dir);
    mG.scale(9.8);
}


void TestRig::constructBase()
{
    Makena::Vec3 p01(-20.0, -2.0,  20.0);
    Makena::Vec3 p02(-20.0,  2.0,  20.0);
    Makena::Vec3 p03( 20.0,  2.0,  20.0);
    Makena::Vec3 p04( 20.0, -2.0,  20.0);
    Makena::Vec3 p05(-20.0, -2.0, -20.0);
    Makena::Vec3 p06(-20.0,  2.0, -20.0);
    Makena::Vec3 p07( 20.0,  2.0, -20.0);
    Makena::Vec3 p08( 20.0, -2.0, -20.0);

    Makena::Manifold m01;

    m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
    Makena::Manifold::Martialled mm01 = m01.exportData();

    Makena::Vec3 com01(0.0, -2.0, 0.0);

    double mass = 1.0;
    double massInv = 1.0/mass;
    double w = 2.0;
    double h = 10.0;
    double d = 2.0;

    Makena::Mat3x3 inertia01( (h*h + d*d) * massInv / 12.0, 0.0, 0.0,
                              0.0, (w*w + d*d) * massInv / 12.0, 0.0,
                              0.0, 0.0, (h*h + w*w) * massInv / 12.0 );

    Makena::Mat3x3 mat01( 1.0, 0.0, 0.0,
                          0.0, 1.0, 0.0,
                          0.0, 0.0, 1.0 );
    double det01 = mat01.det();
    mat01.scale(1.0/det01);
    mBase = new Makena::ConvexRigidBody(0);
    mBase->setObjectAttributes(mm01, mass, inertia01, 10.0, 10.0, true);
    Makena::Quaternion q01(1.0, 0.0, 0.0, 0.0);
    Makena::Vec3 zero(0.0, 0.0, 0.0);
    mBase->setGeomConfig(com01, q01, zero, zero);

    mContM.registerRigidBody(mBase);

}


void TestRig::constructFixed01()
{
    Makena::Vec3 p01(-30.0, -8.0,   5.0);
    Makena::Vec3 p02(-30.0,  8.0,   2.0);
    Makena::Vec3 p03( 30.0,  8.0,   2.0);
    Makena::Vec3 p04( 30.0, -8.0,   5.0);
    Makena::Vec3 p05(-30.0, -8.0,  -5.0);
    Makena::Vec3 p06(-30.0,  8.0,  -5.0);
    Makena::Vec3 p07( 30.0,  8.0,  -5.0);
    Makena::Vec3 p08( 30.0, -8.0,  -5.0);

    Makena::Manifold m01;

    m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
    Makena::Manifold::Martialled mm01 = m01.exportData();

    Makena::Vec3 com01(0.0, 3.0, -21.0);

    double mass = 2.0;
    double massInv = 1.0/mass;
    double w = 2.0;
    double h = 10.0;
    double d = 2.0;

    Makena::Mat3x3 inertia01( (h*h + d*d) * massInv / 12.0, 0.0, 0.0,
                              0.0, (w*w + d*d) * massInv / 12.0, 0.0,
                              0.0, 0.0, (h*h + w*w) * massInv / 12.0 );

    Makena::Mat3x3 mat01( 1.0, 0.0, 0.0,
                          0.0, 1.0, 0.0,
                          0.0, 0.0, 1.0 );
    double det01 = mat01.det();
    mat01.scale(1.0/det01);
    mFixed01 = new Makena::ConvexRigidBody(1);
    mFixed01->setObjectAttributes(mm01, mass, inertia01, 10.0, 10.0, true);
    Makena::Quaternion q01(1.0, 0.0, 0.0, 0.0);
    Makena::Vec3 zero(0.0, 0.0, 0.0);
    mFixed01->setGeomConfig(com01, q01, zero, zero);

    mContM.registerRigidBody(mFixed01);

}


void TestRig::constructFixed02()
{
    Makena::Vec3 p01(-5.0, -8.0,   30.0);
    Makena::Vec3 p02(-2.0,  8.0,   30.0);
    Makena::Vec3 p03( 5.0,  8.0,   30.0);
    Makena::Vec3 p04( 5.0, -8.0,   30.0);
    Makena::Vec3 p05(-5.0, -8.0,  -30.0);
    Makena::Vec3 p06(-2.0,  8.0,  -30.0);
    Makena::Vec3 p07( 5.0,  8.0,  -30.0);
    Makena::Vec3 p08( 5.0, -8.0,  -30.0);

    Makena::Manifold m01;

    m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
    Makena::Manifold::Martialled mm01 = m01.exportData();

    Makena::Vec3 com01(21.0, 3.0, 0.0);

    double mass = 1.0;
    double massInv = 1.0/mass;
    double w = 2.0;
    double h = 10.0;
    double d = 2.0;

    Makena::Mat3x3 inertia01( (h*h + d*d) * massInv / 12.0, 0.0, 0.0,
                              0.0, (w*w + d*d) * massInv / 12.0, 0.0,
                              0.0, 0.0, (h*h + w*w) * massInv / 12.0 );

    Makena::Mat3x3 mat01( 1.0, 0.0, 0.0,
                          0.0, 1.0, 0.0,
                          0.0, 0.0, 1.0 );
    double det01 = mat01.det();
    mat01.scale(1.0/det01);
    mFixed02 = new Makena::ConvexRigidBody(2);
    mFixed02->setObjectAttributes(mm01, mass, inertia01, 10.0, 10.0, true);
    Makena::Quaternion q01(1.0, 0.0, 0.0, 0.0);
    Makena::Vec3 zero(0.0, 0.0, 0.0);
    mFixed02->setGeomConfig(com01, q01, zero, zero);

    mContM.registerRigidBody(mFixed02);

}


static double rand100()
{
    return ((rand()%10000)+1)/100.0;
}

static Makena::Vec3 randVec3D100()
{
    return Makena::Vec3(rand100(), rand100(), rand100());
}



void TestRig::init(){

    constructBase();
    constructFixed01();
    constructFixed02();

    Makena::Vec3 zero(0.0, 0.0, 0.0);

    Makena::Vec3 p01(-1.0, -1.0,  1.0);
    Makena::Vec3 p02(-1.0,  1.0,  1.0);
    Makena::Vec3 p03( 1.0,  1.0,  1.0);
    Makena::Vec3 p04( 1.0, -1.0,  1.0);
    Makena::Vec3 p05(-1.0, -1.0, -1.0);
    Makena::Vec3 p06(-1.0,  1.0, -1.0);
    Makena::Vec3 p07( 1.0,  1.0, -1.0);
    Makena::Vec3 p08( 1.0, -1.0, -1.0);

    Makena::Manifold m01;

    m01.constructCuboid(p01, p02, p03, p04, p05, p06, p07, p08);
    Makena::Manifold::Martialled mm01 = m01.exportData();

    double mass = 1.0;
    double w = 2.0;
    double h = 2.0;
    double d = 2.0;

    Makena::Mat3x3 inertia01( (h*h + d*d) * mass / 12.0, 0.0, 0.0,
                              0.0, (w*w + d*d) * mass / 12.0, 0.0,
                              0.0, 0.0, (h*h + w*w) * mass / 12.0 );

    Makena::Mat3x3 mat01( 1.0, 0.0, 0.0,
                          0.0, 1.0, 0.0,
                          0.0, 0.0, 1.0 );

    double det01 = mat01.det();
    mat01.scale(1.0/det01);

    std::vector<Makena::Vec3> points3D;
    for (long i = 0; i < 10; i++) {
         auto v = randVec3D100();
         v.scale(1.0/50.0);
         points3D.push_back(v);
    }
    Makena::Manifold m02;
    enum Makena::predicate pred;
    m02.findConvexHull(points3D, pred, 1.0e-08);
    Makena::Manifold::Martialled mm02 = m02.exportData();

    mBody01 = new Makena::ConvexRigidBody(101);
    mBody01->setObjectAttributes(mm01, mass, inertia01, 2.0, 2.0, false);
    Makena::Vec3 com01(0.0,  5.0, 0.3);
    Makena::Quaternion q01(cos(M_PI*0.2), 0.0, sin(M_PI*0.2), 0.0);
//    Makena::Quaternion q01(1.0, 0.0, 0.0, 0.0);
    mBody01->setGeomConfig(com01, q01, zero, zero);
    mConsM.registerConvexRigidBody(mBody01);
    mContM.registerRigidBody(mBody01);

    mBody02 = new Makena::ConvexRigidBody(102);
    mBody02->setObjectAttributes(mm01, mass, inertia01, 2.0, 2.0, false);
    Makena::Vec3 com02(0.0,  8.0, 0.2);
//    Makena::Quaternion q02(1.0, 0.0, 0.0, 0.0);
    Makena::Quaternion q02(cos(M_PI*0.3), 0.0, sin(M_PI*0.3), 0.0);
    mBody02->setGeomConfig(com02, q02, zero, zero);
    mConsM.registerConvexRigidBody(mBody02);
    mContM.registerRigidBody(mBody02);

    mBody03 = new Makena::ConvexRigidBody(103);
    mBody03->setObjectAttributes(mm01, mass, inertia01, 2.0, 2.0, false);
    Makena::Vec3 com03(0.0, 11.0, 0.0);
//    Makena::Quaternion q03(1.0, 0.0, 0.0, 0.0);
    Makena::Quaternion q03(cos(M_PI*0.1), 0.0, sin(M_PI*0.1), 0.0);
    mBody03->setGeomConfig(com03, q03, zero, zero);
    mConsM.registerConvexRigidBody(mBody03);
    mContM.registerRigidBody(mBody03);


    mBody04 = new Makena::ConvexRigidBody(104);
    mBody04->setObjectAttributes(mm01, mass, inertia01, 2.0, 2.0, false);
    Makena::Vec3 com04(0.0, 14.0, 0.1);
//    Makena::Quaternion q04(1.0, 0.0, 0.0, 0.0);
    Makena::Quaternion q04(cos(M_PI*0.7), 0.0, sin(M_PI*0.7), 0.0);
    mBody04->setGeomConfig(com04, q04, zero, zero);
    mConsM.registerConvexRigidBody(mBody04);
    mContM.registerRigidBody(mBody04);

    mBody05 = new Makena::ConvexRigidBody(105);
    mBody05->setObjectAttributes(mm01, mass, inertia01, 2.0, 2.0, false);
    Makena::Vec3 com05(0.0, 17.0, 0.2);
    Makena::Quaternion q05(1.0, 0.0, 0.0, 0.0);
    mBody05->setGeomConfig(com05, q05, zero, zero);
    mConsM.registerConvexRigidBody(mBody05);
    mContM.registerRigidBody(mBody05);

    mContM.initializeAABB();

}

void TestRig::term()
{

    if (mBody05!=nullptr) {
        mContM.unregisterRigidBody(mBody05);
        mConsM.unregisterConvexRigidBody(mBody05);
        delete mBody05;
    }
    if (mBody04!=nullptr) {
        mContM.unregisterRigidBody(mBody04);
        mConsM.unregisterConvexRigidBody(mBody04);
        delete mBody04;
    }
    if (mBody03!=nullptr) {
        mContM.unregisterRigidBody(mBody03);
        mConsM.unregisterConvexRigidBody(mBody03);
        delete mBody03;
    }
    if (mBody02!=nullptr) {
        mContM.unregisterRigidBody(mBody02);
        mConsM.unregisterConvexRigidBody(mBody02);
        delete mBody02;
    }
    if (mBody01!=nullptr) {
        mContM.unregisterRigidBody(mBody01);
        mConsM.unregisterConvexRigidBody(mBody01);
        delete mBody01;
    }
    if (mBase!=nullptr) {
        mContM.unregisterRigidBody(mBase);
        delete mBase;
    }
    if (mFixed02!=nullptr) {
        mContM.unregisterRigidBody(mFixed02);
        delete mFixed02;
    }
    if (mFixed01!=nullptr) {
        mContM.unregisterRigidBody(mFixed01);
        delete mFixed01;
    }
}


void TestRig::step(const double deltaT)
{
    if (!mGo) {return;}

    if (mBody01!=nullptr) {
        mBody01->resetExternalForceTorque();
        mBody01->addExternalForce(mG);
    }

    if (mBody02!=nullptr) {
        mBody02->resetExternalForceTorque();
        mBody02->addExternalForce(mG);
    }

    if (mBody03!=nullptr) {
        mBody03->resetExternalForceTorque();
        mBody03->addExternalForce(mG);
    }

    if (mBody04!=nullptr) {
        mBody04->resetExternalForceTorque();
        mBody04->addExternalForce(mG);
    }

    if (mBody05!=nullptr) {
        mBody05->resetExternalForceTorque();
        mBody05->addExternalForce(mG);
    }

    mContM.lockForNextStep(deltaT);
    mContM.updateActiveContacts(mConsM);

    mConsM.initializeStep(deltaT);

    mContM.registerActiveConstraints(mConsM); 

    mConsM.update();

    if (mContM.discoverNewContacts(true)) {

        mContM.registerNewConstraints(mConsM); 

        mConsM.update();
    }

    mConsM.commitUpdate();

    mConsM.terminateStep();
    mContM.unlockAfterStep(mConsM);

}


void TestRig::render(GLBufferManager& bm)
{
    std::vector<Makena::Vec3> verticesFaces;
    std::vector<Makena::Vec3> colorsFaces;
    std::vector<float>        alphasFaces;
    std::vector<Makena::Vec3> normalsFaces;

    std::vector<Makena::Vec3> verticesEdges;
    std::vector<Makena::Vec3> colorsEdges;

    std::vector<Makena::Vec3> verticesPoints;
    std::vector<Makena::Vec3> colorsPoints;


    Makena::Vec3  red(1.0, 0.3, 0.3);
    Makena::Vec3  blue(0.3, 0.3, 1.0);
    Makena::Vec3  gray(0.5, 0.5, 0.5);
    Makena::Vec3  white(1.0, 1.0, 1.0);
    Makena::Vec3  yellow(1.0, 1.0, 0.0);

    verticesEdges.emplace_back(0.0, 0.0, 0.0);
    verticesEdges.push_back(mG);
    colorsEdges.emplace_back(1.0, 1.0, 1.0);
    colorsEdges.emplace_back(1.0, 1.0, 1.0);

    bm.addLines(verticesEdges, colorsEdges);


    if (mBase!=nullptr) {
        Makena::Manifold& m01 = mBase->ConvexHull();
        m01.makeOpenGLVerticesColorsNormalsForTriangles(
            blue, 0.8,
            verticesFaces,
            colorsFaces,
            alphasFaces,
            normalsFaces,
            false,
            mBase->Qmat(),
            mBase->CoM()
        );
    }

    if (mFixed01!=nullptr) {
        Makena::Manifold& m01 = mFixed01->ConvexHull();
        m01.makeOpenGLVerticesColorsNormalsForTriangles(
            blue, 0.8,
            verticesFaces,
            colorsFaces,
            alphasFaces,
            normalsFaces,
            false,
            mFixed01->Qmat(),
            mFixed01->CoM()
        );
    }

    if (mFixed02!=nullptr) {
        Makena::Manifold& m01 = mFixed02->ConvexHull();
        m01.makeOpenGLVerticesColorsNormalsForTriangles(
            blue, 0.8,
            verticesFaces,
            colorsFaces,
            alphasFaces,
            normalsFaces,
            false,
            mFixed02->Qmat(),
            mFixed02->CoM()
        );
    }

    if (mBody01!=nullptr) {
        Makena::Manifold& m01 = mBody01->ConvexHull();
        m01.makeOpenGLVerticesColorsNormalsForTriangles(
            red, 0.8,
            verticesFaces,
            colorsFaces,
            alphasFaces,
            normalsFaces,
            false,
            mBody01->Qmat(),
            mBody01->CoM()
        );
    }

    if (mBody02!=nullptr) {
        Makena::Manifold& m02 = mBody02->ConvexHull();
        m02.makeOpenGLVerticesColorsNormalsForTriangles(
            yellow, 0.8,
            verticesFaces,
            colorsFaces,
            alphasFaces,
            normalsFaces,
            false,
            mBody02->Qmat(),
            mBody02->CoM()
        );
    }

    if (mBody03!=nullptr) {
        Makena::Manifold& m03 = mBody03->ConvexHull();
        m03.makeOpenGLVerticesColorsNormalsForTriangles(
            red, 0.8,
            verticesFaces,
            colorsFaces,
            alphasFaces,
            normalsFaces,
            false,
            mBody03->Qmat(),
            mBody03->CoM()
        );
    }

    if (mBody04!=nullptr) {
        Makena::Manifold& m04 = mBody04->ConvexHull();
        m04.makeOpenGLVerticesColorsNormalsForTriangles(
            yellow, 0.8,
            verticesFaces,
            colorsFaces,
            alphasFaces,
            normalsFaces,
            false,
            mBody04->Qmat(),
            mBody04->CoM()
        );
    }

    if (mBody05!=nullptr) {
        Makena::Manifold& m05 = mBody05->ConvexHull();
        m05.makeOpenGLVerticesColorsNormalsForTriangles(
            red, 0.8,
            verticesFaces,
            colorsFaces,
            alphasFaces,
            normalsFaces,
            false,
            mBody05->Qmat(),
            mBody05->CoM()
        );
    }
    bm.addTriangles(verticesFaces, colorsFaces, alphasFaces, normalsFaces);
}

/***************************************************************************/
/**                   @brief TestRig END                                   */
/***************************************************************************/


static long maxCapacityPoints = 5000;


int main( int argc, char* argv[])
{

    if (argc != 1) {
        std::cerr << "Usage: test_visualizer_joints\n";
        exit(1);
    }

    GlfwManager          glfw   (1024, 768);
    TwManagerAppSpecific tw     (glfw);
    ShaderManager        shaders(tw);
    GLBufferManager      buffers(shaders);
    GCSAxesLines         gcsAxes(5.0, 50.0);
    TestRig              rig;

    if (!glfw.initGLFW()) {
        return 1;
    }

    tw.init();

    glfw.configGLFW();
    glfw.setCallback(&rig);
    tw.setCallback(&rig);

    buffers.genBuffers(maxCapacityPoints);
    shaders.loadShaders();

    rig.init();

    clock_t tPrev, tCur;
    struct timeval tv;
    double microSecondsPrev;
    double microSecondsCur;
    gettimeofday(&tv, NULL);
    microSecondsCur = tv.tv_sec*1000000.0 + tv.tv_usec;
    
    long index = 0;
    do {
        microSecondsPrev = microSecondsCur;
        gettimeofday(&tv, NULL);
        microSecondsCur = tv.tv_sec*1000000.0 + tv.tv_usec;
        
        double deltaT = (microSecondsCur - microSecondsPrev)/1000000.0;
//std::cerr << "deltaT: " << deltaT << " [secs]\n";
//std::cerr << "freq: " << 1.0/deltaT << "[Hz]\n";
        tw.updateVers();

        shaders.setViewGeometry();

        buffers.resetOneStep();

        rig.step(deltaT*2.0);

        buffers.addLines(gcsAxes.mVertices, gcsAxes.mColors);

        rig.render(buffers);

        buffers.draw();

        tw.draw();

        glfw.update();

        
    } while ( glfw.shouldExit() );

    rig.term();

    tw.terminate();
    shaders.unloadShaders();
    buffers.releaseBuffers();
    glfw.terminate();

    return 0;
}

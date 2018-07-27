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
#include "manifold.hpp"
#include "convex_rigid_body.hpp"
#include "binary_dilation.hpp"
#include "gjk_origin_finder.hpp"


class GlfwManager;
class TwManager;

/** @file src_interactive_tests/test_visualizer_gjk_origin_finder.cpp
 *
 *  @brief visual test rig for Makena::GJKOriginFinder
 *
 *  Usage: test_visualizer_gjk_origin_finder <test input file1> 
 *                                            <test input file2>
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


    glm::quat      mBody1ori;
    double         mBody1x;
    double         mBody1y;
    double         mBody1z;

    glm::quat      mBody2ori;
    double         mBody2x;
    double         mBody2y;
    double         mBody2z;

};

TwManagerAppSpecific::TwManagerAppSpecific(GlfwManager&glfw):
    TwManager(glfw),
    mBody1x(0.0),
    mBody1y(0.0),
    mBody1z(0.0),
    mBody2x(0.0),
    mBody2y(0.0),
    mBody2z(0.0){;}



TwManagerAppSpecific::~TwManagerAppSpecific() {;}

void TwManagerAppSpecific::init()
{

    initCommon();

    TwAddVarRW(
        mTwMainBar,
        "Body 1 Orientation",
        TW_TYPE_QUAT4F,
        &mBody1ori,
        " showval=true opened=true "
    );

    TwAddVarRW(
        mTwMainBar,
        "Body 1 X Pos",
        TW_TYPE_DOUBLE,
        &mBody1x,
        " label= 'Body 1 X Pos' min=-20.0 max=20.0 step=0.01 keyIncr=a "
        "keyDecr=A help='Body 1 X Pos' "
    );

    TwAddVarRW(
        mTwMainBar,
        "Body 1 Y Pos",
        TW_TYPE_DOUBLE,
        &mBody1y,
        " label= 'Body 1 Y Pos' min=-20.0 max=20.0 step=0.01 keyIncr=b "
        "keyDecr=B help='Body 1 Y Pos' "
    );

    TwAddVarRW(
        mTwMainBar,
        "Body 1 Z Pos",
        TW_TYPE_DOUBLE,
        &mBody1z,
        " label= 'Body 1 Z Pos' min=-20.0 max=20.0 step=0.01 keyIncr=c "
        "keyDecr=C help='Body 1 Z Pos' "
    );

    TwAddVarRW(
        mTwMainBar,
        "Body 2 Orientation",
        TW_TYPE_QUAT4F,
        &mBody2ori,
        " showval=true opened=true "
    );

    TwAddVarRW(
        mTwMainBar,
        "Body 2 X Pos",
        TW_TYPE_DOUBLE,
        &mBody2x,
        " label= 'Body 2 X Pos' min=-20.0 max=20.0 step=0.01 keyIncr=a "
        "keyDecr=A help='Body 2 X Pos' "
    );

    TwAddVarRW(
        mTwMainBar,
        "Body 2 Y Pos",
        TW_TYPE_DOUBLE,
        &mBody2y,
        " label= 'Body 2 Y Pos' min=-20.0 max=20.0 step=0.01 keyIncr=b "
        "keyDecr=B help='Body 2 Y Pos' "
    );

    TwAddVarRW(
        mTwMainBar,
        "Body 2 Z Pos",
        TW_TYPE_DOUBLE,
        &mBody2z,
        " label= 'Body 2 Z Pos' min=-20.0 max=20.0 step=0.01 keyIncr=c "
        "keyDecr=C help='Body 2 Z Pos' "
    );

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

    TestRig(const long&   maxNumPivots,
            const long&   maxNumCycles,
            const double& epsilonZero   );

    ~TestRig();

    void loadContentsFromFiles(char* filename1, char* filename2);

    void stepMode();

    void immediateMode();

    void GlfwStep();

    void GlfwKeyCallback(char ch);

    void GlfwFindClosestP();

    void TwUpdate(TwManager& tw);

    void render(GLBufferManager& buf);

    std::string strStepType(enum Makena::GJKOriginFinder::stepType t);

  private:
    std::vector<Makena::Vec3> parsePointsFromText(char* filename);

    Makena::ConvexRigidBody               mBody1;
    Makena::ConvexRigidBody               mBody2;
    std::vector<Makena::BDVertex>         mQ;
    Makena::GJKOriginFinder               mFinder;
  
    bool                                  mStepMode;
    bool                                  mResult;

};


TestRig::TestRig(
    const long&   maxNumPivots,
    const long&   maxNumCycles,
    const double& epsilonZero   ):
    mBody1(1),
    mBody2(2),
    mFinder(mBody1, 
            mBody2,
            mQ,
            maxNumPivots,
            maxNumCycles,
            epsilonZero,
            true,
            1.0,
            std::cerr),
    mStepMode(false),
    mResult(false)
{
    Makena::Manifold::Martialled CH;
    Makena::Mat3x3 I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    mBody1.setObjectAttributes(CH, 1.0, I, 0.5, 0.5, false);
    mBody2.setObjectAttributes(CH, 1.0, I, 0.5, 0.5, false);
    mFinder.findOriginSequencerInit(true);
    mFinder.setLogLevel(Makena::GJKOriginFinder::ALL);
}


TestRig::~TestRig(){;}


void TestRig::loadContentsFromFiles(char* filename1, char* filename2)
{
    std::vector<Makena::Vec3> points1 =parsePointsFromText(filename1);
    std::vector<Makena::Vec3> points2 =parsePointsFromText(filename2);
    enum Makena::predicate pred1, pred2;

    mBody1.ConvexHull().findConvexHull(points1, pred1);
    mBody2.ConvexHull().findConvexHull(points2, pred2);
}


void TestRig::stepMode(){ mStepMode = true; }


void TestRig::immediateMode(){;}


void TestRig::GlfwStep()
{
    if (mStepMode) {
        std::cerr<<"One Step BEFORE\n";
        std::cerr<<"    Step No: "  <<mFinder.mST_stepNo               <<"\n";
        std::cerr<<"    Step Type: "<<strStepType(mFinder.mST_stepType)<<"\n";

        mFinder.findOriginSequencerOneStep();

        if (mFinder.mST_stepType == Makena::GJKOriginFinder::ST_RETURN_TYPE1) {
            mResult = true;
        }
        else if( ( mFinder.mST_stepType == 
                   Makena::GJKOriginFinder::ST_RETURN_TYPE2 ) ||
                 ( mFinder.mST_stepType == 
                   Makena::GJKOriginFinder::ST_RETURN_TYPE2 )   ) {
            mResult = false;
        }
        std::cerr << "One Step AFTER\n";
    }
}


void TestRig::GlfwKeyCallback(char ch)
{
    if (ch == 's') {
        // Snapshot.
        std::cerr << "Body 1\n";
        auto vpair1 = mBody1.ConvexHull().vertices();
        for (auto vit = vpair1.first; vit != vpair1.second; vit++) {
            auto p = (*vit)->pGCS(mBody1.Qmat(), mBody1.CoM());
            std::cerr << p << "\n";
        }
        std::cerr << "Body 2\n";
        auto vpair2 = mBody2.ConvexHull().vertices();
        for (auto vit = vpair2.first; vit != vpair2.second; vit++) {
            auto p = (*vit)->pGCS(mBody2.Qmat(), mBody2.CoM());
            std::cerr << p << "\n";
        }
    }
}


void TestRig::GlfwFindClosestP(){
    if (!mStepMode) {
        mResult = mFinder.findOrigin(true);
    }
}


void TestRig::TwUpdate(TwManager& tw)
{
    auto& tw2 = dynamic_cast<TwManagerAppSpecific&>(tw);

    Makena::Vec3 zero(0.0, 0.0, 0.0);
    Makena::Quaternion q1( tw2.mBody1ori.w,
                           tw2.mBody1ori.x,
                           tw2.mBody1ori.y,
                           tw2.mBody1ori.z  );
    q1.normalize();
    Makena::Vec3 CoM1(tw2.mBody1x, tw2.mBody1y, tw2.mBody1z);
    mBody1.setGeomConfig(CoM1, q1, zero, zero);
    mBody1.updateGeomConfig(0.0);

    Makena::Quaternion q2( tw2.mBody2ori.w,
                           tw2.mBody2ori.x,
                           tw2.mBody2ori.y,
                           tw2.mBody2ori.z  );
    q2.normalize();
    Makena::Vec3 CoM2(tw2.mBody2x, tw2.mBody2y, tw2.mBody2z);
    mBody2.setGeomConfig(CoM2, q2, zero, zero);
    mBody2.updateGeomConfig(0.0);
}

std::vector<Makena::Vec3> TestRig::parsePointsFromText(char* filename)
{
    long numPoints = 0;
    std::vector<Makena::Vec3> points;

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
            if (numPoints++ > 100) {
                break;
            }
        }
    }
    gs.close();
    return points;
}


std::string TestRig::strStepType(enum Makena::GJKOriginFinder::stepType t)
{
    switch (t) {
      case Makena::GJKOriginFinder::ST_INIT:
        return "ST_INIT";
        break;

      case Makena::GJKOriginFinder::ST_BEFORE_loop:
        return "ST_BEFORE_loop";
        break;

      case Makena::GJKOriginFinder::ST_BEGIN_loop:
        return "ST_BEGIN_loop";
        break;

      case Makena::GJKOriginFinder::ST_AFTER_findClosestPointToOrigin:
        return "ST_AFTER_findClosestPointToOrigin";
        break;

      case Makena::GJKOriginFinder::ST_AFTER_removeUnsupportingVertices:
        return "ST_AFTER_removeUnsupportingVertices";
        break;

      case Makena::GJKOriginFinder::
                               ST_AFTER_findExtremePointsAlongSupportDirection:
        return " ST_AFTER_findExtremePointsAlongSupportDirection";
        break;

      case Makena::GJKOriginFinder::ST_AFTER_updateSimplexTowardOrigin:
        return "ST_AFTER_updateSimplexTowardOrigin";
        break;

      case Makena::GJKOriginFinder::ST_AFTER_loop:
        return "ST_AFTER_loop";
        break;

      case Makena::GJKOriginFinder::ST_RETURN_TYPE1:
        return "ST_RETURN_TYPE1";
        break;

      case Makena::GJKOriginFinder::ST_RETURN_TYPE2:
        return "ST_RETURN_TYPE2";
        break;

      case Makena::GJKOriginFinder::ST_RETURN_TYPE3:
        return "ST_RETURN_TYPE3";
        break;

      case Makena::GJKOriginFinder::ST_RETURN_TYPE4:
        return "ST_RETURN_TYPE4";
        break;

      case Makena::GJKOriginFinder::ST_RETURN_TYPE5:
        return "ST_RETURN_TYPE5";
        break;

      case Makena::GJKOriginFinder::ST_TERMINATED:
        return "ST_TERMINATED";
        break;
    }
    return "";
}


void TestRig::render(GLBufferManager& buf)
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

    mFinder.renderObj1Surface(
          red,  0.8, verticesFaces, colorsFaces, alphasFaces, normalsFaces);
    mFinder.renderObj2Surface(
          blue, 0.8, verticesFaces, colorsFaces, alphasFaces, normalsFaces);
    mFinder.renderBDManifoldSurface(
          gray, 0.5, verticesFaces, colorsFaces, alphasFaces, normalsFaces);
    buf.addTriangles(verticesFaces, colorsFaces, alphasFaces, normalsFaces);
                     
    mFinder.renderObj1Lines         (red,  verticesEdges, colorsEdges);
    mFinder.renderObj2Lines         (blue, verticesEdges, colorsEdges);
    mFinder.renderBDManifoldLines(gray, verticesEdges, colorsEdges);
    mFinder.renderSupportVectorLine(yellow, verticesEdges, colorsEdges);
    buf.addLines(verticesEdges, colorsEdges);

    mFinder.renderActiveFeature1Points(red,  verticesPoints, colorsPoints);
    mFinder.renderActiveFeature2Points(blue, verticesPoints, colorsPoints);
    mFinder.renderActiveFeatureBDManifoldPoints(
                                       gray, verticesPoints, colorsPoints);

    mFinder.renderClosestPoint
                      (mResult?yellow:white, verticesPoints, colorsPoints);

    buf.addPoints(verticesPoints, colorsPoints);


}

/***************************************************************************/
/**                   @brief TestRig END                                   */
/***************************************************************************/


static long maxCapacityPoints = 5000;


int main( int argc, char* argv[])
{

    if (argc != 3&&argc != 4) {
        std::cerr << "Usage: test_visualizer_gjk_origin_finder "
                     "[test input file A] [test input file B] [STEP]\n";
        exit(1);
    }

    GlfwManager          glfw   (1024, 768);
    TwManagerAppSpecific tw     (glfw);
    ShaderManager        shaders(tw);
    GLBufferManager      buffers(shaders);
    GCSAxesLines         gcsAxes(5.0, 50.0);
    TestRig              rig    (100, 8, 0.0);

    rig.loadContentsFromFiles(argv[1], argv[2]);

    if (argc == 4) {
        if (strcmp(argv[3], "STEP")==0) {
            rig.stepMode();
        }
    }

    if (!glfw.initGLFW()) {
        return 1;
    }

    tw.init();

    glfw.configGLFW();
    glfw.setCallback(&rig);
    tw.setCallback(&rig);

    buffers.genBuffers(maxCapacityPoints);
    shaders.loadShaders();

    do {

        tw.updateVers();

        shaders.setViewGeometry();

        buffers.resetOneStep();

        rig.GlfwFindClosestP();

        // Add triangles, lines, and points here.

        buffers.addLines(gcsAxes.mVertices, gcsAxes.mColors);
        rig.render(buffers);

        buffers.draw();

        tw.draw();

        glfw.update();

    } while ( glfw.shouldExit() );

    tw.terminate();
    shaders.unloadShaders();
    buffers.releaseBuffers();
    glfw.terminate();

    return 0;
}

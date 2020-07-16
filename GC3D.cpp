//#include "trackball.h"
//#include "GC3D.h"
//#include "cell.h"
//#include "lattice.h"
//#include <iostream>
//#include <iomanip> // setprecision
////#include "glm/glm/glm.hpp"
//#include <sstream>
//#ifdef WIN32
//#include "GL/glut.h"
//#endif
//#ifdef __linux__
//#include <GL/glut.h>
//#endif
//#ifdef __APPLE__
//#include <OpenGL/gl.h>
//#include <OpenGL/glu.h>
//#include <GLUT/glut.h>
//#endif
//using namespace std;
//#define NB_TYPES_CELLS 6
//
//vector<B_cell*>* currentBC_list;
//vector<T_cell*>* currentTC_list;
//vector<FDC*>* currentFDC_list;
//vector<Plasma_cell*>* currentPC_list;
//vector<Memory_cell*>* currentMEM_list;
//vector<Stromal_cell*>* currentSTROMA_list;
//lattice* currentSpace;
//double currentTime;
//bool is_moving = 0;
//int begining_x, begining_y;
//int W = 400, H = 600;
//float curent_quaternion[4];
//float last_quaternion[4];
//bool reModel = 1;
//float scalefactor = 1.0;
//
//
//
//void review(int w, int h)
//{
//    glViewport(0, 0, w, h);
//    W = w;
//    H = h;
//}
//void recalcModelView(void)
//{
//    GLfloat m[4][4];
//    glPopMatrix();
//    glPushMatrix();
//    build_rotmatrix(m, curent_quaternion);
//    glMultMatrixf(&m[0][0]);
//    if (scalefactor == 1.0) {
//        glDisable(GL_NORMALIZE);
//    } else {
//        glEnable(GL_NORMALIZE);
//    }
//    glScalef(scalefactor, scalefactor, scalefactor);
//    reModel = 0;
//}
//
//void mouse(int button, int state, int x, int y)
//{
//    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && not(pause1)) {
//        glutIdleFunc(NULL);
//        is_moving = 1;
//        begining_x = x;
//        begining_y = y;
//    }
//    if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
//        is_moving = 0;
//    }
//    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && pause1)
//    {
//        GLdouble posX, posY, posZ,posXs, posYs, posZs,posXe, posYe, posZe;
//        glutIdleFunc(NULL);
//        GLfloat window_width = glutGet(GLUT_WINDOW_WIDTH);
//        GLfloat window_height = glutGet(GLUT_WINDOW_HEIGHT);
//        GLbyte color[4];
//        GLfloat depth;
//        GLuint index;
//        int d = window_width < window_height ? window_width : window_height;            /* minimum dimension    */
//        int xl = ( window_width - d ) / 2;
//        int yb = ( window_height - d ) / 2;
//        
//        GLint viewport[4];
//        glGetIntegerv(GL_VIEWPORT, viewport);
//        
//        GLdouble modelview[16];
//        glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
//        GLdouble projection[16];
//        glGetDoublev(GL_PROJECTION_MATRIX, projection);
//        double winX = x;
//        double winY = viewport[3] - (double)y;
//        glLoadIdentity();
//        glReadPixels(x, winY, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, color);
//        glReadPixels(x, winY , 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
//        glReadPixels(x,winY, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &index);
//        printf("Clicked on pixel %d, %d, color %02hhx%02hhx%02hhx%02hhx, depth %f, stencil index %u\n",
//               x, y, color[0], color[1], color[2], color[3], depth, index);
//        gluUnProject( x, winY, depth, modelview, projection, viewport, &posX, &posY, &posZ);
//        printf("Coordinates in object space: %f, %f, %f\n",
//               &posX, &posY, &posZ);
//        //        is_moving = 1;
//        
//    }
//    if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
//        is_moving = 0;
//    }
//    if(button == 3){
//        scalefactor *= 1.08;
//        reModel = 1;
//        glutPostRedisplay();
//    }
//    if(button == 4){
//        scalefactor *= 0.92;
//        reModel = 1;
//        glutPostRedisplay();
//    }
//}
//
//void motion(int x, int y)
//{
//    
//
//    if (is_moving) {
//        trackball(last_quaternion,
//                  (2.0 * begining_x - W) / W,
//                  (H - 2.0 * begining_y) / H,
//                  (2.0 * x - W) / W,
//                  (H - 2.0 * y) / H
//                  );
//        begining_x = x;
//        begining_y = y;
//        add_quats(last_quaternion, curent_quaternion, curent_quaternion);
//        reModel = 1;
//        glutPostRedisplay();
//    }
//}
//
//
//
//
//
//
//
//
//
//
//
//
//void drawBitmapText(string s,float x,float y,float z)
//{
//    glRasterPos3f(x, y,z);
//    for(unsigned int i = 0; i < s.size(); ++i)
//        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, (int) s[i]);
//}
//
//void drawBitmapText2D(string s,float x,float y)
//{
//    glRasterPos2f(x, y);
//    for(unsigned int i = 0; i < s.size(); ++i)
//        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, (int) s[i]);
//}
//
//void keyboard(unsigned char key, int x, int y){
//    switch(key) {
//        case 'w':
//            scalefactor *= 1.08;
//            reModel = 1;             // needs to recompute the model (object)
//            glutPostRedisplay(); // needs
//            break;
//        case 's':
//            scalefactor *= 0.92;
//            reModel = 1;
//            glutPostRedisplay();
//            break;
//        default:
//            break;
//    }
//    glutPostRedisplay();
//}
//unsigned int PickBuffer[512];
//void init(void);
//void initGC3D(int argc, char** argv)
//{
//    
//    glutInit(&argc, argv);
//    
//    // Says the type of buffer (memory) used
//    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE| GLUT_STENCIL);
//    glutInitWindowSize(600, 600);
//    glutInitWindowPosition(0,0);
//    // initialize the trackball (calculations for how much you turn the object based on mouse movements)
//    trackball(curent_quaternion, 0.0, 0.0, 0.0, 0.0);
//    glutCreateWindow("GC3D");
//    glRenderMode(GL_SELECT);
//    glutDisplayFunc(display);     // you give the function to refresh the picture
//    glutReshapeFunc(review);   // you give the function to refresh when the size of screen is changed
//    glutMouseFunc(mouse);
////    glutMouseWheelFunc(mouse); //Elena: 06-05-2019 glutMouseWheelFunc(mouse) returns segfault
//    glutMotionFunc(motion);
//    glutKeyboardFunc(keyboard);
//    glutCreateMenu(controlLights);
//    glutAddMenuEntry("Pause", 3);
//    glutAddMenuEntry("Quit", 4);
//    glutAttachMenu(GLUT_RIGHT_BUTTON);
//    glEnable(GL_CULL_FACE);
////  glSelectBuffer( 512, PickBuffer );
//    
//    
//    //#ifdef __APPLE__
//    //    #ifdef USE_IL_LIBRARY
//    //    ilInit();
//    //    iluInit();
//    //    ilutRenderer(ILUT_OPENGL);
//    //    #endif
//    //#endif
//    
//    init();
//    
//}
//
//
//void init(void)
//{
//    currentBC_list = NULL;
//    currentTC_list = NULL;
//    currentFDC_list = NULL;
//    currentPC_list = NULL;
//    currentMEM_list = NULL;
//    currentSTROMA_list = NULL;
//    glClearStencil(0); // this is the default value
//    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
//    glEnable(GL_STENCIL_TEST);
//    glClearColor(0.0, 0.0, 0.0, 1.0);
//    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
//    glShadeModel(GL_SMOOTH);
//    glEnable(GL_LIGHTING);
//    glEnable(GL_LIGHT0);
//    GLfloat lightOnePosition[] = {0, 0, 100, 0.0};
//    glLightfv(GL_LIGHT0, GL_POSITION, lightOnePosition);
//    GLfloat lightOneColor[] = {1.0, 1.0, 1.0, 1.0};
//    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightOneColor);
//#define zoom 1
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    gluPerspective(60.0, 1.0, 1, 500);
//    glMatrixMode(GL_MODELVIEW);
//    gluLookAt(0, 0, 65* 2.0,0, 0, 0,0.0, 1.0, 0.);
//    glEnable(GL_COLOR_MATERIAL);
//    glColorMaterial(GL_FRONT, GL_DIFFUSE);
//    glEnable(GL_DEPTH_TEST);
//}
//
//void nextToDisplay(vector<B_cell*>* BC_list, vector<T_cell*> *TC_list,
//                   vector<FDC*>* FDC_list, vector<Plasma_cell*>* PC_list, vector<Memory_cell*>* MEM_list, vector<Stromal_cell*>* STROMA_list, lattice* s, double t){
//    currentBC_list = BC_list;
//    currentTC_list = TC_list;
//    currentFDC_list = FDC_list;
//    currentPC_list = PC_list;
//    currentMEM_list = MEM_list;
//    currentSTROMA_list = STROMA_list;
//    currentSpace = s;
//    currentTime = t;
//}
//
//
//void clearDisplay(){
//    currentBC_list = NULL;
//    currentTC_list = NULL;
//    currentFDC_list = NULL;
//    currentPC_list = NULL;
//    currentMEM_list = NULL;
//    currentSTROMA_list = NULL;
//    currentSpace = NULL;
//    currentTime = 0;
//}
//
//
//void display(){
//#define XWidth currentSpace->X
//#define YWidth currentSpace->Y
//#define ZWidth currentSpace->Z
//#define sizeSphere 0.4
//    if (reModel)
//        recalcModelView();
//    int w = glutGet(GLUT_WINDOW_WIDTH);
//    int h = glutGet(GLUT_WINDOW_HEIGHT);
//    
//    //    glMatrixMode( GL_PROJECTION );
//    //    glLoadIdentity();
//    //
//    //    glm::vec3 cameraPos = glm::vec3(1, 1, 1);
//    //    glm::vec3 cameraTarget = glm::vec3(0,0,0);
//    //    glm::vec3 cameraDirection = glm::normalize(cameraPos - cameraTarget);
//    //    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);
//    //    glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
//    //    glm::vec3 cameraRight = glm::normalize(glm::cross(up, cameraDirection));
//    //    glm::vec3 cameraUp = glm::cross(cameraDirection, cameraRight);
//
//    //
//    //    if (currentTime>0.1 && currentTime<0.2){
//    //
//    //    gluLookAt(cameraPos.x,cameraPos.y,cameraPos.z,cameraPos.x+cameraFront.x,cameraPos.y+cameraFront.y,cameraPos.z+cameraFront.z,cameraUp.x,cameraUp.y,cameraUp.z);
//    //        glutPostRedisplay();
//    //        reModel=1;
//    //    }
//    //    else {
//    //    gluLookAt(10,cameraPos.y,cameraPos.z,cameraPos.x+cameraFront.x,cameraPos.y+cameraFront.y,cameraPos.z+cameraFront.z,cameraUp.x,cameraUp.y,cameraUp.z);
//    //        glutPostRedisplay();
//    //        reModel=1;
//    //    }
//    
//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
//    glClearStencil(0);
//    glEnable(GL_STENCIL_TEST);
//    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
//    glMatrixMode(GL_MODELVIEW);
//    glLineWidth(3.0);
//    glColor4f(0.45, 0.0, 0.45, 1.0);
//    glColor4f(1, 1, 1, 1.0);
//    for(int ix = 0; ix <=1; ++ix){
//        for(int iy = 0; iy <=1; ++iy){
//            for(int iz = 0; iz <=1; ++iz){
//                glBegin(GL_LINES);
//                glVertex3f(XWidth * (ix-0.5) , YWidth * (iy-0.5) , -ZWidth * (iz-0.5) );
//                glVertex3f(XWidth * (ix-0.5) , YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
//                glVertex3f(XWidth * (ix-0.5) , YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
//                glVertex3f(XWidth * (ix-0.5) , - YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
//                glVertex3f(XWidth * (ix-0.5) , - YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
//                glVertex3f(XWidth * (ix-0.5) , YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
//                glVertex3f(XWidth * (ix-0.5) , YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
//                glVertex3f(- XWidth * (ix-0.5) , YWidth * (iy-0.5) , ZWidth * (iz-0.5) );
//                glEnd();
//            }
//        }
//    }
//    
//    // Axes
//    glColor4f(0.3, 0.3, 0.85, 1.0);
//    for(int i = 0; i < XWidth ; ++i){
//        glBegin(GL_LINES);
//        glVertex3f(i -0.1 - 0.5*XWidth, 0 , 0 );
//        glVertex3f(i + 0.1 - 0.5*XWidth, 0 , 0 );
//        glEnd();
//        glBegin(GL_LINES);
//        glVertex3f(0,i - 0.1 - 0.5*YWidth, 0  );
//        glVertex3f(0,i + 0.1 - 0.5*YWidth, 0 );
//        glEnd();
//        glBegin(GL_LINES);
//        glVertex3f(0,0,i - 0.1  - 0.5*ZWidth);
//        glVertex3f(0,0,i + 0.1 - 0.5*ZWidth);
//        glEnd();
//    }
//    drawBitmapText("X",XWidth / 2.0,0,0);
//    drawBitmapText("Y",0,YWidth / 2.0,0);
//    drawBitmapText("Z",0,0,ZWidth/2.0);
//    auto x_str = std::to_string(currentTime);
//    
//    drawBitmapText("t= "+x_str+" (hrs)",0,0,ZWidth/2.0+20);
//    
//    //Danial: Last update 21-10-18 01:46, checked--> has a problem with glutsSolidSphere
//    bool show_chemo=false;
//    if (show_chemo)
//    {
//        for (int i=0; i<65;i++)
//        {
//            for (int j=0; j<65;j++)
//            {
//                for (int k=0; k<65;k++)
//                {
//                    double cl12 =currentSpace->chemoat(CXCL12, i, j, k)/301.103;
//                    double cl13 =currentSpace->chemoat(CXCL13, i, j, k)/4.51654;
//                    glColor4f(cl12, 0., cl13,1); //red for cl13
//                    glTranslatef(i - XWidth / 2., j -YWidth / 2., k - ZWidth / 2.);
//                    glutSolidSphere(0.1, 5, 5);
//                    glTranslatef(-i + XWidth / 2., -j +YWidth / 2., -k + ZWidth / 2.);
//                }
//            }
//        }
//        
//    }
//    if(currentBC_list){
//        for(unsigned int is = 0; is < currentBC_list->size(); ++is){
//            B_cell* currCB = (*currentBC_list)[is];
//            if(currCB->cell_type == Centroblast)
//                glColor4f(0.2, 0.2, 1.0, 1.0);
//            if(currCB->cell_type == Centrocyte)
//                glColor4f(0.97, 0.97, 0.05, 1.0);
//            glTranslatef(currCB->position.X - XWidth / 2., currCB->position.Y -YWidth / 2., currCB->position.Z - ZWidth / 2.);
//            glutSolidSphere(sizeSphere, 20, 20);
//            glTranslatef(-currCB->position.X + XWidth / 2., -currCB->position.Y +YWidth / 2., -currCB->position.Z + ZWidth / 2.);
//            glStencilFunc(GL_ALWAYS, currCB->ID, -1);
//        }
//    }
//    
//    if(true) if(currentFDC_list){
//        for(unsigned int is = 0; is < currentFDC_list->size(); ++is){
//            
//            FDC* currFDC = (*currentFDC_list)[is];
//            glLoadName( currFDC->ID );
//            
//            glColor4f(1.0, 0.5, 0.0, 1.0);
//            glTranslatef(currFDC->position.X - XWidth / 2., currFDC->position.Y -YWidth / 2., currFDC->position.Z - ZWidth / 2.);
//            glutSolidCube(sizeSphere);         // Philippe I think this way you can draw cubes!!
//            glTranslatef(-currFDC->position.X + XWidth / 2., -currFDC->position.Y +YWidth / 2., -currFDC->position.Z + ZWidth / 2.);
//            
//            // plotting tails of the FDC in another color / size
//            int V = currFDC->occupiedPositions.size();
//            for(int fr = 0; fr < V; ++fr){
//                
//                glColor4f(1.0, 0.5, 0, 1.0);
//                glTranslatef(currFDC->occupiedPositions[fr].X - XWidth / 2., currFDC->occupiedPositions[fr].Y -YWidth / 2., currFDC->occupiedPositions[fr].Z - ZWidth / 2.);
//                glutSolidCube(sizeSphere / 5.0);
//                glTranslatef(-currFDC->occupiedPositions[fr].X  + XWidth / 2., -currFDC->occupiedPositions[fr].Y +YWidth / 2., -currFDC->occupiedPositions[fr].Z  + ZWidth / 2.);
//                glStencilFunc(GL_ALWAYS, currFDC->ID, -1);
//                
//            }
//            glStencilFunc(GL_ALWAYS, currFDC->ID, -1);
//            
//        }
//    }
//    
//    
//    if(currentTC_list){
//        for(unsigned int is = 0; is < currentTC_list->size(); ++is){
//            T_cell* currTC = (*currentTC_list)[is];
//            // cout << "T cell" << is << " at pos " << currTC->position.X << "," << currTC->position.Y << "," << currTC->position.Z << endl;
//            
//            glLoadName( currTC->ID );
//            if (currTC->ID==617)
//            {
//                float x,y,z;
//                x=currTC->position.X;
//                y=currTC->position.Y;
//                z=currTC->position.Z;
//                
//                glColor4d(0.0, 1.0, 0.1, 1.0);
//                glTranslatef(currTC->position.X - XWidth / 2., currTC->position.Y -YWidth / 2., currTC->position.Z - ZWidth / 2.);
//                glutWireSphere(sizeSphere, 20, 20);
//                glTranslatef(-currTC->position.X + XWidth / 2., -currTC->position.Y +YWidth / 2., -currTC->position.Z + ZWidth / 2.);
//            }
//            else
//            {
//                glEnable(GL_BLEND);
//                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//                glColor4f(0.0, 1.0, 0.1, 0.4);
//                glTranslatef(currTC->position.X - XWidth / 2., currTC->position.Y -YWidth / 2., currTC->position.Z - ZWidth / 2.);
//                glutSolidSphere(sizeSphere, 20, 20);
//                glTranslatef(-currTC->position.X + XWidth / 2., -currTC->position.Y +YWidth / 2., -currTC->position.Z + ZWidth / 2.);
//            }
//            glStencilFunc(GL_ALWAYS, currTC->ID, 1);
//        }
//    }
//    
//    //    if(currentSTROMA_list){
//    //        for(unsigned int is = 0; is < currentSTROMA_list->size(); ++is){
//    //            Stromal_cell* currSTR = (*currentSTROMA_list)[is];
//    //            glColor4f(1.0, 0.5, 0.0, 5.0);
//    //            glTranslatef(currSTR->position.X - XWidth / 2., currSTR->position.Y -YWidth / 2., currSTR->position.Z - ZWidth / 2.);
//    //            glutSolidSphere(sizeSphere, 20, 20);
//    //            glTranslatef(-currSTR->position.X + XWidth / 2., -currSTR->position.Y +YWidth / 2., -currSTR->position.Z + ZWidth / 2.);
//    //        }
//    //    }
//    
//    if(currentPC_list){
//        for(unsigned int is = 0; is < currentSTROMA_list->size(); ++is){
//            Plasma_cell* Plasma = (*currentPC_list)[is];
//            if (not(Plasma==NULL)){
//                if (not(Plasma->cell_state==Plasma_Out)){
//                    glColor4f(1.0, 1., 1., 5.0);
//                    glTranslatef(Plasma->position.X - XWidth / 2., Plasma->position.Y -YWidth / 2., Plasma->position.Z - ZWidth / 2.);
//                    glutSolidSphere(sizeSphere, 20, 20);
//                    glTranslatef(-Plasma->position.X + XWidth / 2., -Plasma->position.Y +YWidth / 2., -Plasma->position.Z + ZWidth / 2.);
//                }
//            }
//        }
//    }
//    glutSwapBuffers();  // in the double buffers mode (see init)
//#ifdef USE_IL_LIBRARY
//    vector< unsigned char > buf( w * h * 3 );
//    glPixelStorei( GL_PACK_ALIGNMENT, 1 );
//    glReadPixels( 0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &buf[0] );
//    static int cpt = 1000;
//    cpt ++;
//    stringstream newFile;
//    
//#endif
//}
//
//static void idle(void)
//{
//    static float time = 0.0;
//    static float jump = 0.0;
//
//    time = glutGet(GLUT_ELAPSED_TIME) / 500.0;
//
//    jump = 4.0 * fabs(sin(time)*0.5);
////    if (!lightMoving) {
////        lightAngle += 0.03;
////    }
//    glutPostRedisplay();
//}
//
//GLboolean lightZeroSwitch = GL_TRUE, lightOneSwitch = GL_TRUE;
//
//void
//controlLights(int value)
//{
//    switch (value) {
//        case 3:
//            pause1= 1-pause1;
//            //            if (glIsEnabled(GL_MULTISAMPLE_ARB)) {
//            //                glDisable(GL_MULTISAMPLE_ARB);
//            //            } else {
//            //                glEnable(GL_MULTISAMPLE_ARB);
//            //            }
//            break;
//        case 4:
//            glutFullScreen();
//            break;
//        case 5:
//            exit(0);
//            break;
//    }
//    glutPostRedisplay();
//}
//
//
//void processHits (GLint hits, GLuint buffer[])
//{
//    unsigned int i, j;
//    GLuint names, *ptr;
//    
//    printf ("hits = %d\n", hits);
//    ptr = (GLuint *) buffer;
//    for (i = 0; i < hits; i++) { /*  for each hit  */
//        names = *ptr;
//        printf (" number of names for hit = %d\n", names); ptr++;
//        printf("  z1 is %g;", (float) *ptr/0x7fffffff); ptr++;
//        printf(" z2 is %g\n", (float) *ptr/0x7fffffff); ptr++;
//        printf ("   the name is ");
//        for (j = 0; j < names; j++) {     /*  for each name */
//            printf ("%d ", *ptr); ptr++;
//        }
//        printf ("\n");
//    }
//}
//
//#define BUFSIZE 512
//void selectObjects(void)
//{
//    GLuint selectBuf[BUFSIZE];
//    GLint hits;
//    
//    glSelectBuffer (BUFSIZE, selectBuf);
//    (void) glRenderMode (GL_SELECT);
//    
//    glInitNames();
//    glPushName(0);
//    
//    glPushMatrix ();
//    glMatrixMode (GL_PROJECTION);
//    glLoadIdentity ();
//    glOrtho (0.0, 5.0, 0.0, 5.0, 0.0, 10.0);
//    glMatrixMode (GL_MODELVIEW);
//    glLoadIdentity ();
//    glLoadName(1);
//    glPopMatrix ();
//    glFlush ();
//    
//    hits = glRenderMode (GL_RENDER);
//    processHits (hits, selectBuf);
//}
//
//

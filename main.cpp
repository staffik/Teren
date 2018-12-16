#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>

#include <common/shader.hpp>

#include <fstream>
#include "dirent.h"
#include <regex>

using namespace std;
using namespace glm;

#define READ_ONLY_NEEDED 1
#define LOGGING_ENABLED 1

const int TILE_SIDE = 1201;
const float EARTH_RADIUS = 2;

#ifdef LOGGING_ENABLED
  #define LOG(msg) cout<<msg<<endl;
#else
  #define LOG(msg)
#endif

#define XY2I(i,j) ((i)*TILE_SIDE+(j))

ostream& operator<<(ostream &os, const glm::vec3 &v) {
  os<<"("<<v.x<<" "<<v.y<<" "<<v.z<<")";
  return os;
}

ostream& operator<<(ostream &os, const glm::mat4 m) {
  for(int i=0; i<4; ++i, cout<<endl) for(int j=0; j<4; ++j) cout<<m[i][j]<<" ";
}

float POSITIONS[2*TILE_SIDE*TILE_SIDE];
float NET[8];
GLuint POSITIONS_VBO;
GLuint NET_VBO;

glm::mat4 zero_mat() {
	glm::mat4 m;
	for(int i=0; i<4; ++i) for(int j=0; j<4; ++j) m[i][j]=0;
	return m;
}

glm::mat4 id_mat() {
	glm::mat4 res(zero_mat());
	for(int i=0; i<4; ++i) res[i][i]=1;
	return res;
}

glm::vec3 spherical_to_cartesian(float a, float b) {
  glm::vec3 res;
  res.x = cos(a) * cos(b);
  res.y = sin(b);
  res.z = sin(a) * cos(b);
  return res;
}

glm::vec3 local_to_global(glm::vec3 v, float lat, float lon) {
  glm::mat4 m = id_mat();
  m = glm::rotate(m, glm::radians(-lat), glm::vec3(0.0f,1.0f,0.0f));
  m = glm::rotate(m, glm::radians(-90+lon), glm::vec3(0.0f,0.0f,1.0f));
  glm::vec3 o = m * glm::vec4(v, 1.0);
  return o;
}

struct coord {
  int longitude;
  int latitude;

  coord() = default;

  coord(int x, int y) {
    longitude = y;
    latitude = x;
  }

  pair<int,int> to_pair() const {
    return make_pair(longitude, latitude);
  }

  bool operator<(const coord &C) const {
    return to_pair() < C.to_pair();
  }
};

struct tile {
  coord c;
  vector<float> h;
  GLuint VAO;
  GLuint VBO;
  int elements_no;

  void init() {
    assert(h.size() == TILE_SIDE*TILE_SIDE);
    
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, POSITIONS_VBO);
    glVertexAttribPointer(0, 1, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), (GLvoid*)0);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), (GLvoid*)(sizeof(GLfloat)));
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
   
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(h[0])*h.size(), h.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(2);

    glBindVertexArray(0);
  }

	void init_net() {
    elements_no = 8;
    for(int i=0; i<elements_no; ++i) h.push_back(-1);
 		glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, NET_VBO);
    glVertexAttribPointer(0, 1, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), (GLvoid*)0);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), (GLvoid*)(sizeof(GLfloat)));
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);

    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(h[0])*h.size(), h.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(2);

    glBindVertexArray(0);
  }
};

tile net_tile;

struct lod {
  int level;
  GLuint EBO;
  vector<unsigned int> indices;
  int elements_no;

  void init(int lvl) {
    assert(lvl>0 && lvl<10);
    level = lvl;
    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);

    int skip = (1<<lvl);
    int hs = (skip >> 1);
    int n = (TILE_SIDE >> lvl);

    for (int i=0; i<TILE_SIDE+skip; i+=skip) {
      for (int j=0; j<TILE_SIDE+skip; j+=skip) {
        int ii = std::min(i, TILE_SIDE-1);
        int jj = std::min(j, TILE_SIDE-1);

        int a = std::min(ii+hs, TILE_SIDE-1);
        int b = std::min(jj+hs, TILE_SIDE-1);
        int c = std::max(jj-hs, 0);
        indices.push_back(XY2I(ii,jj));
        indices.push_back(XY2I(a,b));
        indices.push_back(XY2I(a,c));

        a = std::max(ii-hs, 0);
        b = std::min(jj+hs, TILE_SIDE-1);
        c = std::max(jj-hs, 0);
        indices.push_back(XY2I(ii,jj));
        indices.push_back(XY2I(a,b));
        indices.push_back(XY2I(a,c));

        a = std::max(jj-hs, 0);
        b = std::min(ii+hs, TILE_SIDE-1);
        c = std::max(ii-hs, 0);
        indices.push_back(XY2I(ii,jj));
        indices.push_back(XY2I(b,a));
        indices.push_back(XY2I(c,a));

        a = std::min(jj+hs, TILE_SIDE-1);
        b = std::min(ii+hs, TILE_SIDE-1);
        c = std::max(ii-hs, 0);
        indices.push_back(XY2I(ii,jj));
        indices.push_back(XY2I(b,a));
        indices.push_back(XY2I(c,a));
      }
    }
    elements_no = indices.size();

    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices[0])*elements_no, indices.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  }

	void init_net() {
    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);

    for(int i=0; i<4; ++i)
      indices.push_back(i);
    elements_no = indices.size();

    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices[0])*elements_no, indices.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  }
};

lod net_lod;

struct camera {
  enum class VIEW_MODE { MAP=0, GLOBE };

  float SPEED_3D = 0.00005;
  float SPEED_MAP = 0.001;
  float MAX_ABS_PITCH=89.0;
  float MOUSE_SENSIVITY=0.1;
  float MIN_RAD = 0.1;
  double MAX_FOV=60;

  glm::vec3 INIT_POS_2D   = glm::vec3(0.0, 0.0, 3.0);
  glm::vec3 INIT_FRONT_2D = glm::vec3(0.0, 1.0, 0.0);
  glm::vec3 INIT_UP_2D    = glm::vec3(0.0, 0.0, 1.0);

  float INIT_RAD_3D   = EARTH_RADIUS*1.01; // dist from earth's center
  glm::vec3 INIT_POS_3D   = glm::vec3(0.0, INIT_RAD_3D, 0.0);
  glm::vec3 INIT_FRONT_3D = glm::vec3(1.0, 0.0, 0.0);
  glm::vec3 INIT_UP_3D    = glm::vec3(0.0, 1.0, 0.0);
  float INIT_ALPHA_3D = glm::radians(0.0f); // latitude
  float INIT_BETA_3D  = glm::radians(90.0f); // longitude

  VIEW_MODE view_mode = VIEW_MODE::MAP;

  glm::vec3 pos = INIT_POS_2D;
  glm::vec3 front = INIT_FRONT_2D;
  glm::vec3 up = INIT_UP_2D;
  float alpha;
  float beta;
  float rad;

  float speed = SPEED_MAP;
  float yaw_value   = 0;
  float pitch_value =  0.0;

  int SCR_WIDTH=1024;
  int SCR_HEIGHT=768;
  bool mouse_started = true;
  float fov_value   =  45.0;
  double previous_time;
  float previous_x=0.0;
  float previous_y=0.0;

  glm::mat4 projection_mat;
  glm::mat4 view_mat;

  void init() {
    if (view_mode == VIEW_MODE::MAP) {
      pos = INIT_POS_2D;
      front = INIT_FRONT_2D;
      up = INIT_UP_2D;
			speed = SPEED_MAP;
    } else {
      alpha = INIT_ALPHA_3D;
      beta = INIT_BETA_3D;
      rad = INIT_RAD_3D;
      pos = INIT_POS_3D;
      front = INIT_FRONT_3D;
      up = INIT_UP_3D;
			speed = SPEED_3D;
    }
  }

  void update_front_vec() {
    glm::vec3 v = spherical_to_cartesian(glm::radians(yaw_value), glm::radians(pitch_value));
    v = local_to_global(v, glm::degrees(alpha), glm::degrees(beta));
    front = glm::normalize(v);
  }

  void adjust(float a, float b) {
    a = glm::radians(a);
    b = glm::radians(b);
    alpha = a; beta = b;
    up = spherical_to_cartesian(a, b);
    pos = up * rad;
    update_front_vec();
  }

  void recalculate_geometry() {
    projection_mat = glm::perspective(glm::radians(fov_value), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.001f, 10000.0f);
    view_mat = glm::lookAt(pos, pos + front, up);
  }
} Camera;

struct state {
  string path;
  DIR *dir = NULL;
  GLFWwindow *window = NULL;
  GLuint programID_2D = 0;
  GLuint programID_3D = 0;
  GLuint u_transform_2D;
  GLuint u_transform_3D;
  GLuint u_model;
  GLuint u_R;

  int min_x = -180;
  int max_x = 180;
  int min_y = -90;
  int max_y = 90;

  int beg_x = -180;
  int end_x = 180;
  int beg_y = -90;
  int end_y = 90;

  float H = 180.0;
  float W = 360.0;
  float h = 2.0/H;
  float w = 2.0/W;
  float SCALE_HEIGHT = 1;
  float SCALE_X = 1.0;
  float SCALE_Y = 1.0;


  bool auto_lod = false;
  int current_lod = 1;
  float previous_time;
  long long triangles_counter = 0;
  int frames_counter = 0;

  int MAX_FPS = 20;

  map<coord, tile> Tiles;
  map<int, lod> LoDs;

  GLuint uniform_transform() {
    return Camera.view_mode == camera::VIEW_MODE::MAP ? u_transform_2D : u_transform_3D;
  }

  float middle_x() {
    return (min_x+max_x)*0.5;
  }

  float middle_y() {
    return (min_y+max_y)*0.5;
  }
  
} State;

glm::mat4 get_model_mat(int x, int y) {
  float w = State.w, h = State.h;
  if (Camera.view_mode == camera::VIEW_MODE::GLOBE) {
    w = 2.0/State.W;
    h = 2.0/State.H;
  }
  glm::vec3 translation((x-State.beg_x)*w, (y-State.beg_y)*h, 0.0);
  glm::mat4 M = id_mat();
  M = glm::translate(M, glm::vec3(-1.0, -1.0, 0.0));
  M = glm::translate(M, translation);
  M = glm::scale(M, glm::vec3(w/1200.0, h/1200.0, 1.0));

  return M;
}

void update_ranges() {
  if (Camera.view_mode == camera::VIEW_MODE::MAP) {
    State.beg_x = State.min_x;
    State.end_x = State.max_x;
    State.beg_y = State.min_y;
    State.end_y = State.max_y;
  } else {
    State.beg_x = -180;
    State.end_x = 180;
    State.beg_y = -90;
    State.end_y = 90;
  }

  float scale_x = cos((State.beg_y+State.end_y)*M_PI/360.0);
  assert(scale_x > 0);

  State.H = State.end_y - State.beg_y + 1;
  State.W = State.end_x - State.beg_x + 1;
  State.h = std::min(2.0 / State.H, 2.0 / (State.W * scale_x));
  State.w = State.h * scale_x;
}

void set_3d_uniforms() {
  glUniform1f(State.u_R, EARTH_RADIUS);
}

void draw() {
  Camera.recalculate_geometry();
  lod &lo = State.LoDs[State.current_lod];
  glm::mat4 view_mat(id_mat());
  glm::mat4 projection_mat(id_mat());

  if (Camera.view_mode == camera::VIEW_MODE::GLOBE) {
    view_mat = Camera.view_mat;
    projection_mat = Camera.projection_mat;
  } else {
    view_mat = glm::lookAt(glm::vec3(Camera.pos.x, Camera.pos.y, 3), glm::vec3(Camera.pos.x, Camera.pos.y, 0.0), glm::vec3(0.0, 1.0, 0.0));  
    projection_mat = Camera.projection_mat;
  }

	//for (int x=State.beg_x; x<=State.end_x; ++x) {
    //for (int y=State.beg_y; y<=State.end_y; ++y) {
  for (int x=State.min_x; x<=State.max_x; ++x) {
    for (int y=State.min_y; y<=State.max_y; ++y) {
           glm::mat4 model_mat = get_model_mat(x, y);
      glm::mat4 transform_mat = projection_mat * view_mat;
      if (Camera.view_mode == camera::VIEW_MODE::GLOBE) {
        glUniformMatrix4fv(State.u_model, 1, GL_FALSE, &model_mat[0][0]);
        set_3d_uniforms();
      }
      else {
        transform_mat *= model_mat;
      }
      glUniformMatrix4fv(State.uniform_transform(), 1, GL_FALSE, &transform_mat[0][0]);

			auto it = State.Tiles.find(coord(x,y));
      if (it == State.Tiles.end()) {
				glBindVertexArray(net_tile.VAO);
				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, net_lod.EBO);
        glDrawElements(GL_LINE_LOOP, net_tile.elements_no, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);
      	State.triangles_counter += net_tile.elements_no;
			} else {
				tile &t = it->second;
				glBindVertexArray(t.VAO);
				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, lo.EBO);
				glDrawElements(GL_TRIANGLES, lo.elements_no, GL_UNSIGNED_INT, 0);
				glBindVertexArray(0);
				State.triangles_counter += lo.elements_no / 3;
			}
    }
  }
  State.frames_counter++;
}

void set_lod(int L) {
   if (L == 0) {
     State.auto_lod = true;
   } else {
     State.auto_lod = false;
     State.current_lod = L;
   }
}

void change_position(glm::vec3 delta_pos) {
  assert(Camera.view_mode == camera::VIEW_MODE::MAP);
  Camera.pos += glm::vec3(delta_pos.x, delta_pos.y, 0);
}

void change_position_3d(glm::vec2 delta_pos) {
  assert(Camera.view_mode == camera::VIEW_MODE::GLOBE);
  Camera.alpha += delta_pos.x;
  Camera.beta += delta_pos.y;
  Camera.beta = std::min(Camera.beta, (float)M_PI_2);
  Camera.beta = std::max(Camera.beta, (float)-M_PI_2);

  Camera.up = spherical_to_cartesian(Camera.alpha, Camera.beta);
  Camera.pos = Camera.rad * Camera.up;

  Camera.update_front_vec();
}

void get_front_right_vecs(glm::vec2 &front, glm::vec2 &right) {
  glm::vec2 v(1,0);
  v = glm::rotate(v, glm::radians(-90+Camera.yaw_value));
  front = v;
  v = glm::rotate(v, glm::radians(90.0f));
  right = v;
}

void keyboard_handling() {
  float time_delta = glfwGetTime() - Camera.previous_time;
  float camera_delta = Camera.speed * time_delta;

  if (Camera.view_mode == camera::VIEW_MODE::MAP)
    camera_delta *= 1.0 / sqrt(Camera.fov_value);

  if (Camera.view_mode == camera::VIEW_MODE::GLOBE) {
    glm::vec2 front, right;
    get_front_right_vecs(front, right);
    front *= camera_delta;
    right *= camera_delta;

    if (glfwGetKey(State.window, GLFW_KEY_PAGE_UP) == GLFW_PRESS) {
        Camera.rad += camera_delta * EARTH_RADIUS;
        Camera.pos = Camera.up * Camera.rad;
    }

    if (glfwGetKey(State.window, GLFW_KEY_PAGE_DOWN) == GLFW_PRESS) {
        Camera.rad -= camera_delta * EARTH_RADIUS;
        Camera.rad = std::max(Camera.rad, Camera.MIN_RAD);
        Camera.pos = Camera.up * Camera.rad;
    }

    if (glfwGetKey(State.window, GLFW_KEY_W) == GLFW_PRESS) {
        glm::vec2 delta_pos = front;
        change_position_3d(delta_pos);
    }

    if (glfwGetKey(State.window, GLFW_KEY_A) == GLFW_PRESS) {
        glm::vec2 delta_pos = -right;
        change_position_3d(delta_pos);
    }

    if (glfwGetKey(State.window, GLFW_KEY_S) == GLFW_PRESS) {
        glm::vec2 delta_pos = -front;
        change_position_3d(delta_pos);
    }

    if (glfwGetKey(State.window, GLFW_KEY_D) == GLFW_PRESS) {
        glm::vec2 delta_pos = right;
        change_position_3d(delta_pos);
    }
  }
  else {
    if (glfwGetKey(State.window, GLFW_KEY_W) == GLFW_PRESS) {
        glm::vec3 delta_pos = camera_delta * Camera.front;
        change_position(delta_pos);
    }

    if (glfwGetKey(State.window, GLFW_KEY_A) == GLFW_PRESS) {
        glm::vec3 delta_pos = -glm::normalize(glm::cross(Camera.front, Camera.up)) * camera_delta;
        change_position(delta_pos);
    }

    if (glfwGetKey(State.window, GLFW_KEY_S) == GLFW_PRESS) {
        glm::vec3 delta_pos = -camera_delta * Camera.front;
        change_position(delta_pos);
    }

    if (glfwGetKey(State.window, GLFW_KEY_D) == GLFW_PRESS) {
        glm::vec3 delta_pos = glm::normalize(glm::cross(Camera.front, Camera.up)) * camera_delta;
        change_position(delta_pos);
    }
  }

  if (glfwGetKey(State.window, GLFW_KEY_0) == GLFW_PRESS) set_lod(0);
  if (glfwGetKey(State.window, GLFW_KEY_1) == GLFW_PRESS) set_lod(1);
  if (glfwGetKey(State.window, GLFW_KEY_2) == GLFW_PRESS) set_lod(2);
  if (glfwGetKey(State.window, GLFW_KEY_3) == GLFW_PRESS) set_lod(3);
  if (glfwGetKey(State.window, GLFW_KEY_4) == GLFW_PRESS) set_lod(4);
  if (glfwGetKey(State.window, GLFW_KEY_5) == GLFW_PRESS) set_lod(5);
  if (glfwGetKey(State.window, GLFW_KEY_6) == GLFW_PRESS) set_lod(6);
  if (glfwGetKey(State.window, GLFW_KEY_7) == GLFW_PRESS) set_lod(7);
  if (glfwGetKey(State.window, GLFW_KEY_8) == GLFW_PRESS) set_lod(8);
  if (glfwGetKey(State.window, GLFW_KEY_9) == GLFW_PRESS) set_lod(9);

  static int view_mode_toggle = GLFW_RELEASE;

  int view_mode_key = glfwGetKey(State.window, GLFW_KEY_TAB);
  if (view_mode_key == GLFW_RELEASE && view_mode_toggle == GLFW_PRESS) {
      if (Camera.view_mode == camera::VIEW_MODE::MAP) {
          Camera.view_mode = camera::VIEW_MODE::GLOBE;
          update_ranges();
          Camera.init();
          Camera.adjust(State.middle_x(), State.middle_y());
          glUseProgram(State.programID_3D);
      } else {
          Camera.view_mode = camera::VIEW_MODE::MAP;
          update_ranges();
          Camera.init();
          glUseProgram(State.programID_2D);
      }
  }
  view_mode_toggle = view_mode_key;
}

void change_viewsize(GLFWwindow *window, int width, int height) {
  Camera.SCR_WIDTH=width;
  Camera.SCR_HEIGHT=height;
  glViewport(0, 0, width, height);
}

void moving_mouse_handling(GLFWwindow* window, double pos_x, double pos_y) {
  if (Camera.view_mode == camera::VIEW_MODE::MAP)
    return;

  if (Camera.mouse_started) {
    Camera.previous_x = pos_x;
    Camera.previous_y = pos_y;
    Camera.mouse_started = false;
  }

  float x_offset = pos_x - Camera.previous_x;
  float y_offset = Camera.previous_y - pos_y;

  x_offset *= Camera.MOUSE_SENSIVITY;
  y_offset *= Camera.MOUSE_SENSIVITY;

  Camera.yaw_value += x_offset;
  Camera.pitch_value += y_offset;

  Camera.previous_x = pos_x;
  Camera.previous_y = pos_y;

  Camera.pitch_value = std::max(-Camera.MAX_ABS_PITCH,
                          std::min(Camera.MAX_ABS_PITCH, Camera.pitch_value));

  Camera.update_front_vec();
}

void scrolling_handling(GLFWwindow* window, double xoffset, double yoffset) {
  Camera.fov_value -= yoffset;
	Camera.fov_value = std::max(1.0, std::min(Camera.MAX_FOV, (double)Camera.fov_value));
}

int gl_stuff_init() {
	if (!glfwInit()) {
		cerr<<"Failed to initialize GLFW."<<endl;
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); 
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	State.window = glfwCreateWindow(1024, 768, "Teren", NULL, NULL);

	if (State.window == NULL) {
		cerr << "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials." << endl;
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(State.window);
	glfwSetFramebufferSizeCallback(State.window, change_viewsize);

	glewExperimental = true;
	if (glewInit() != GLEW_OK) {
		cerr<<"Failed to initialize GLEW."<<endl;
		glfwTerminate();
		return -1;
	}
  
  glfwSetCursorPosCallback(State.window, moving_mouse_handling);
  glfwSetScrollCallback(State.window, scrolling_handling);
	glfwSetInputMode(State.window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwSetInputMode(State.window, GLFW_STICKY_KEYS, GL_TRUE);

	State.programID_2D = LoadShaders( "vs_2d.vertexshader", "FragmentShader.fragmentshader" );
	State.programID_3D = LoadShaders( "vs_3d.vertexshader", "FragmentShader.fragmentshader" );
	glUseProgram(State.programID_2D);

  for(int i=1; i<10; ++i)
    State.LoDs[i].init(i);

  glGenBuffers(1, &POSITIONS_VBO);
  glBindBuffer(GL_ARRAY_BUFFER, POSITIONS_VBO);

  for (int y=TILE_SIDE-1, i=0; y>=0; --y)
    for (int x=0; x<TILE_SIDE; ++x) {
      POSITIONS[i++] = x;
      POSITIONS[i++] = y;
    }

  glBufferData(GL_ARRAY_BUFFER, sizeof(POSITIONS), POSITIONS, GL_STATIC_DRAW);

	glGenBuffers(1, &NET_VBO);
  glBindBuffer(GL_ARRAY_BUFFER, NET_VBO);

	NET[0] = 0;
	NET[1] = 0;
	NET[2] = 0;
	NET[3] = TILE_SIDE;
	NET[4] = TILE_SIDE;
	NET[5] = TILE_SIDE;
	NET[6] = TILE_SIDE;
	NET[7] = 0;

  glBufferData(GL_ARRAY_BUFFER, sizeof(NET), NET, GL_STATIC_DRAW);

  State.u_transform_2D = glGetUniformLocation(State.programID_2D, "u_transform");
  State.u_transform_3D = glGetUniformLocation(State.programID_3D, "u_transform");
  State.u_model = glGetUniformLocation(State.programID_3D, "u_model");
  State.u_R = glGetUniformLocation(State.programID_3D, "u_R");

  glEnable(GL_DEPTH_TEST); 

  return 0;
}

bool out_of_interest(coord c) {
  int y = c.longitude, x = c.latitude;
  return x < State.min_x || x > State.max_x || y < State.min_y || y > State.max_y;
}

int read_tile(string name, tile &t) {
  regex r("([NS])([0-9]{2})([EW])([0-9]{3}).hgt");
  smatch m;
  regex_match(name, m, r);
  if (m.size() != 5)
    return -1;

  int y=stoi(m[2].str().c_str());
  if (m[1]=="S") y=-y;
  int x = stoi(m[4].str().c_str());
  if (m[3]=="W") x=-x;

  t.c.latitude = x;
  t.c.longitude = y;

  if (READ_ONLY_NEEDED && out_of_interest(t.c))
    return -1;
  
  ifstream f(State.path+name, ios::in | ios::binary);
  if (!f.is_open()) {
    cerr << "Error opening file: "<<State.path+name << endl;
    return -1;
  }

  int col_idx = 0;
  short h, last_in_row = 0;
  vector<short> last_in_col(TILE_SIDE, 0);
  while (f.read(reinterpret_cast<char *> (&h), sizeof(short))) {
    h=__builtin_bswap16(h);
    
    if (h == -32768)
      h = ((int)last_in_row + (int)last_in_col[col_idx])/2;
    last_in_row = h;
    last_in_col[col_idx] = h;
    col_idx++; if (col_idx==TILE_SIDE) col_idx=0;
   
    t.h.push_back(h);
  }

  f.close();
  
  return 0;
}

int arg_error(string msg) {
  cerr<<endl;
  cerr<<msg<<endl;
  cerr<<"Example usage:\n"
        "./main.out path/to/hgt/data/ -sz 73 74 -dl 42 43\n"
        "or\n"
        "./main.out path/to/hgt/data/"<<endl;
  cerr<<endl;
  return -1;
}

int read_rng(char* argv[], int from) {
  int mode = 0;
  if (string(argv[from]) == "-sz") mode = 2;
  if (string(argv[from]) == "-dl") mode = 1;
  if (mode == 0) return -1;
  int beg = stoi(argv[from+1]);
  int end = stoi(argv[from+2])-1;
  if (mode==1) { State.min_x=beg; State.max_x=end; }
  if (mode==2) { State.min_y=beg; State.max_y=end; }
  return 0;
}

int parse_options(int argc, char* argv[]) {
  if (argc !=2 && argc!=5 && argc!=8) return arg_error("Invalid number of arguments");
  if ((State.dir = opendir(argv[1])) == NULL) return arg_error("Could not open the directory");
  State.path = string(argv[1])+"/";
  for (int i=0; i<(argc-2)/3; ++i) {
    if (read_rng(argv, 2+i*3) != 0) return arg_error("Could not parse latitude/longitude");
  }
}

int read_data() {
  LOG("READING TILES ...");
  int cnt_readed = 0;
  dirent *ent;
  while ((ent = readdir(State.dir)) != NULL) {
    string name(ent->d_name);
    if (name.size()<4 || name.substr(name.size()-4)!=".hgt")
      continue;

    tile t;
    if (read_tile(name, t) == 0) {
      t.init();
      State.Tiles[t.c] = t;
      LOG("Read tile: "+name);
      cnt_readed++;
    }
  }
  closedir(State.dir);

  if (cnt_readed == 0) {
    cerr << "NO TILES READED";
    return -1;
  }
  LOG("DONE READING "<<cnt_readed<<" TILES!");

  return 0;
}

void clear_line() {
  cout << "\r";
  for (int i=0; i<80; ++i) cout<<" ";
  cout << "\r";
}


void calculate_geometry() {
  auto it = State.Tiles.begin();
  State.min_x = State.max_x = it->first.latitude;
  State.min_y = State.max_y = it->first.longitude;

  for (auto it : State.Tiles) {
    State.min_x = std::min(State.min_x, it.first.latitude);
    State.max_x = std::max(State.max_x, it.first.latitude);
    State.min_y = std::min(State.min_y, it.first.longitude);
    State.max_y = std::max(State.max_y, it.first.longitude);
  }

  update_ranges();
  set_3d_uniforms();
}

void update_state(float current_time) {
  State.previous_time = current_time;
  if (State.auto_lod) {
    if (State.frames_counter <= 10 && State.current_lod+1 < 10)
      State.current_lod++;
    if (State.current_lod > 1 && State.frames_counter >= 60)
      State.current_lod--;
  }
  clear_line();
  cout << State.triangles_counter << "\ttriangles / second\t";
  cout << "|||\t" << State.frames_counter << " FPS\t";
  cout << "|||\tLOD: " << State.current_lod;
  cout.flush();
  State.triangles_counter = State.frames_counter = 0;
}

int main(int argc, char* argv[]) {
  if (parse_options(argc, argv) == -1)
    return EXIT_FAILURE;

  if (gl_stuff_init() == -1)
		return EXIT_FAILURE;

  if (read_data() == -1)
    return EXIT_FAILURE;

  calculate_geometry();
	net_tile.init_net();
	net_lod.init_net();

 	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  double previous_time = State.previous_time = glfwGetTime();
	while (glfwGetKey(State.window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
        !glfwWindowShouldClose(State.window)) {

    double current_time=glfwGetTime();
    if (current_time - State.previous_time > 1.0)
      update_state(current_time);
    
    if (State.MAX_FPS>0 && current_time-previous_time<1.0/State.MAX_FPS) continue;
    previous_time = current_time;
		
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    keyboard_handling();
    draw();
    
    glfwSwapBuffers(State.window);
    glfwPollEvents();
	}

	glfwTerminate();
  cout << endl;

	return EXIT_SUCCESS;
}


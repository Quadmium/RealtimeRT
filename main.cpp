#include <GL/glew.h>
#include <GL/glu.h>
#include <SDL2/SDL.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <chrono>
#include <unordered_map>

#include "vec3.h"

std::chrono::time_point<std::chrono::high_resolution_clock> Now() {
  return std::chrono::high_resolution_clock::now();
}

double ElapsedSeconds(
    std::chrono::time_point<std::chrono::high_resolution_clock> start,
    std::chrono::time_point<std::chrono::high_resolution_clock> end) {
  return static_cast<std::chrono::duration<double>>(end - start).count();
}

const char* glErrorString(GLenum err) {
  switch (err) {
    case GL_NO_ERROR:
      return "GL_NO_ERROR: No error has been recorded";
    case GL_INVALID_ENUM:
      return "GL_INVALID_ENUM: An unacceptable value is specified for an "
             "enumerated argument";
    case GL_INVALID_OPERATION:
      return "GL_INVALID_OPERATION: The specified operation is not allowed in "
             "the current state";
    case GL_INVALID_FRAMEBUFFER_OPERATION:
      return "GL_INVALID_FRAMEBUFFER_OPERATION: The command is trying to "
             "render "
             "to or read from the framebuffer while the currently bound "
             "framebuffer is not framebuffer complete";
    case GL_OUT_OF_MEMORY:
      return "GL_OUT_OF_MEMORY: There is not enough memory left to execute the "
             "command";
  }
  return "Unspecified Error";
}

// Simple check error call
int check_gl_error(const char* call) {
  int err = glGetError();
  if (err != 0) {
    int prog;
    glGetIntegerv(GL_CURRENT_PROGRAM, &prog);
    fprintf(stderr, "%s\nGlError: '%s' CurrentProgram: %i\n\n", call,
            glErrorString(err), prog);
    exit(1);
  }
  return err;
}

// We can wrap the check error in a define to disable the checking
#define GL_DEBUG(m) \
  m;                \
  check_gl_error(#m);

// TODO Figure out how to get this compound statement to work so it's one r
// value statment. Current breaks on void return values. #define GL_DEBUG(m)
//({__auto_type ret = m; check_gl_error(#m); ret;})

// Compiles a shader
static GLint compile_shader(GLenum shader_type, const char* shader_file) {
  // Open our shader file
  FILE* f;
  f = fopen(shader_file, "r");
  if (!f) {
    fprintf(stderr, "Unable to open shader file %s. Aborting.\n", shader_file);
    return -1;
  }
  // Get the shader size
  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  fprintf(stderr, "Compiling %s. Shader size: %li bytes\n", shader_file, fsize);
  char shader_source[fsize + 1];
  fseek(f, 0, SEEK_SET);
  fread(shader_source, fsize, 1, f);
  shader_source[fsize] = '\0';
  fclose(f);
  // Compile the shader
  GLuint shader = GL_DEBUG(glCreateShader(shader_type));
  const char* ss = shader_source;
  glShaderSource(shader, 1, &ss, NULL);
  glCompileShader(shader);
  // Check status
  GLint status;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
  if (status == GL_FALSE) {
    // Output error if there was a problem
    GLint info_log_length;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &info_log_length);
    GLchar info_log[info_log_length];
    glGetShaderInfoLog(shader, info_log_length, NULL, info_log);
    fprintf(stderr, "Compile failure in shader:\n%s\n", info_log);
    return -1;
  }
  return shader;
}

// Verticies for our output triangle
//
// We enable culling as an example below, so  the order these
// are provided actually matter. We need, using the right hand
// rule, for the trangle's face to be towards the screen
//
// Recall that screen space coordinates are laid out like the following:
//  0.0  0.0 -> Center
// -1.0 -1.0 -> Bottom left
//  1.0  1.0 -> Top right
static const GLfloat triangle[][2] = {
    {0.0f, 4.0f},    // Top Middle
    {-4.0f, -4.0f},  // Bottom Left
    {4.0f, -4.0f}    // Bottom Right
};

// Color for our output triangle, corresponding to our verticies
// In this case R, G, B, and an Alpha channel
static const GLfloat colors[][4] = {
    {1.0, 0.0, 0.0, 1.0},  // Red
    {0.0, 1.0, 0.0, 1.0},  // Green
    {0.0, 0.0, 1.0, 1.0},  // Blue
};

int main(int argc, char** argv) {
  // Initialize SDL
  if (SDL_Init(SDL_INIT_EVERYTHING) < 0) {
    fprintf(stderr, "Unable to initialize SDL\n");
    fprintf(stderr, "%s\n", SDL_GetError());
    return 1;
  }
  // Initialize OpenGL context
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 5);
  // Enabled multisampling (a type of antialiasing0
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
  SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
  SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 8);
  // Create the SDL window
  SDL_Window* window = SDL_CreateWindow(
      "Hello Triangle", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 600,
      600, SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL);
  if (window == NULL) {
    fprintf(stderr, "Unable to create SDL Window\n");
    fprintf(stderr, "%s\n", SDL_GetError());
    exit(1);
  }
  // Create the SDL window
  SDL_GLContext gl_context = SDL_GL_CreateContext(window);
  if (gl_context == NULL) {
    fprintf(stderr, "Unable to create OpenGL Context\n");
    fprintf(stderr, "%s\n", SDL_GetError());
    exit(1);
  }
  // Output OpenGL Version
  const unsigned char* version = GL_DEBUG(glGetString(GL_VERSION));
  if (version == NULL) {
    fprintf(stderr, "Unable to determine OpenGL version\n");
    exit(1);
  }
  fprintf(stderr, "OpenGl Version: %s\n", version);
  // Set current window
  SDL_GL_MakeCurrent(window, gl_context);
  // Initialize glew (which loads the OpenGL Function Pointers)
  GLenum glew_status = glewInit();
  if (glew_status != 0) {
    fprintf(stderr, "Error: %s\n", glewGetErrorString(glew_status));
    exit(1);
  }

  // enable depth testing. This doesn't matter for a single triangle but
  // useful if you have a 3d scene or multiple triangles that overlap
  //
  // While you could get away with rendering triangles in the right order,
  // this is difficult to do in complex scenes or scenes where you need
  // per-pixel ordering.
  /*GL_DEBUG(glEnable(GL_DEPTH_TEST));
  GL_DEBUG(glDepthFunc(GL_LESS));
  GL_DEBUG(glDepthRange(0, 1000));
  // Enable alpha blending
  GL_DEBUG(glEnable(GL_BLEND));
  GL_DEBUG(glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA));
  GL_DEBUG(glEnable(GL_MULTISAMPLE));*/
  // We can cull triangles that are wound away from us. Also important for 3d
  // scenes but not this so much. Note that since we have it enabled here
  // the order we provide the verticies in the triangle array matters
  //
  // The rule that if the triangle is wound facing the viewer it will be shown
  //GL_DEBUG(glEnable(GL_CULL_FACE));

  GL_DEBUG(glEnable(GL_BLEND));
  // Compile our shaders
  GLint vertex_shader;
  GLint fragment_shader;
  vertex_shader = compile_shader(GL_VERTEX_SHADER, "vertex.glsl");
  fragment_shader = compile_shader(GL_FRAGMENT_SHADER, "fragment.glsl");
  // If we have a problem quit
  if (vertex_shader < 0 || fragment_shader < 0)
    exit(1);

  // Link the shaders together into a single 'shader program'
  GLuint program;
  program = glCreateProgram();
  glAttachShader(program, vertex_shader);
  glAttachShader(program, fragment_shader);
  glLinkProgram(program);
  // error checking
  GLint status;
  // Get link log if there is a problem
  glGetProgramiv(program, GL_LINK_STATUS, &status);
  if (status == GL_FALSE) {
    GLint info_log_length;
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &info_log_length);

    GLchar info_log[info_log_length];
    glGetProgramInfoLog(program, info_log_length, NULL, info_log);
    fprintf(stderr, "Shader linker failure: %s\n", info_log);
    exit(1);
  }
  // Once the shaders are linked into a program they don't need to be kept
  // around
  glDeleteShader(vertex_shader);
  glDeleteShader(fragment_shader);
  // Handles for various buffers
  // vao is a buffer of buffers, (so will point our verticies and colors)
  // vbo_verticies is for the triangle points and
  // vbo_colors is for the triangle colors
  GLuint vao, vbo_verticies, vbo_colors;

  // Generate a single Vertex Array and bind it  (meaning subsequent calls to do
  // with vertex arrays will refer to it)
  GL_DEBUG(glGenVertexArrays(1, &vao));
  GL_DEBUG(glBindVertexArray(vao));
  // Generate the buffer object for the verticies
  GL_DEBUG(glGenBuffers(1, &vbo_verticies));
  // bind it for the next few calls
  GL_DEBUG(glBindBuffer(GL_ARRAY_BUFFER, vbo_verticies));
  // Upload our triangle data to the vbo
  GL_DEBUG(glBufferData(GL_ARRAY_BUFFER, sizeof(triangle), triangle,
                        GL_STATIC_DRAW));
  // We get the location of the 'in_position' named in the vertex shader
  GLint in_position_loc = GL_DEBUG(glGetAttribLocation(program, "in_position"));
  // Set the location in the vao to this buffer and tell it how to access the
  // data. We have 2 points per vertex hence 2, and sizeof(float) * 2 and the
  // GL_FLOAT
  GL_DEBUG(glVertexAttribPointer(in_position_loc, 2, GL_FLOAT, GL_FALSE,
                                 sizeof(float) * 2, 0));
  // Enable this buffer
  GL_DEBUG(glEnableVertexAttribArray(in_position_loc));
  // Now geneate the vbo for colors
  //GL_DEBUG(glGenBuffers(1, &vbo_colors));
  // Bind it for the next few calls
  //GL_DEBUG(glBindBuffer(GL_ARRAY_BUFFER, vbo_colors));
  // Upload the color data in the same way we did triangles
  //GL_DEBUG(
  //    glBufferData(GL_ARRAY_BUFFER, sizeof(colors), colors, GL_STATIC_DRAW));
  // We get the location of the 'in_color' named in the vertex shader
  //GLint in_color_loc = GL_DEBUG(glGetAttribLocation(program, "in_color"));
  // This time we have RGBA values so set up 4 floats per vertex
  //GL_DEBUG(glVertexAttribPointer(in_color_loc, 4, GL_FLOAT, GL_FALSE,
  //                               sizeof(float) * 4, 0));
  // Enable the vbo
  //GL_DEBUG(glEnableVertexAttribArray(in_color_loc));
  // Now we set to use the shader program we previously compiled
  glUseProgram(program);

  GLint g_seed = glGetUniformLocation(program, "g_seed");
  GLint u_blend = glGetUniformLocation(program, "blend");
  GLint u_scene = glGetUniformLocation(program, "scene");
  glUniform1i(u_scene, 1);
  GLint u_camera_pos = glGetUniformLocation(program, "camera_pos");
  GLint u_camera_forward = glGetUniformLocation(program, "camera_forward");
  GLint u_camera_up = glGetUniformLocation(program, "camera_up");
  GLint u_camera_right = glGetUniformLocation(program, "camera_right");

  vec3 camera_pos = {0, 0, 1};
  vec3 camera_forward = {0, 0, -1};

  auto last_time = Now();
  auto start_time = Now();

  double total_time = 0;
  double last_total_time = 0;

  float cam_angle_up = -3.14159/2;

  SDL_GL_SetSwapInterval(0);
  std::unordered_map<int, bool> keys;

  long long frames_still = 0;

  while(true) {
    auto cur_time = Now();
    double delta_time = ElapsedSeconds(last_time, cur_time);
    last_time = cur_time;

    last_total_time = total_time;
    total_time = ElapsedSeconds(start_time, cur_time);

    if (int(total_time) % 2 != int(last_total_time) % 2) {
      printf("FPS ~ %f\n", 1/delta_time);
    }

    vec3 camera_right = cross(camera_forward, {0, 1, 0});
    vec3 camera_up = cross(camera_right, camera_forward);
    bool moving = false;

    if (keys[SDLK_w]) {
      camera_pos += camera_forward * delta_time;
      moving = true;
    }
    if (keys[SDLK_s]) {
      camera_pos += -1 * camera_forward * delta_time;
      moving = true;
    }
    if (keys[SDLK_a]) {
      camera_pos += -1 * camera_right * delta_time;
      moving = true;
    }
    if (keys[SDLK_d]) {
      camera_pos += 1 * camera_right * delta_time;
      moving = true;
    }
    if (keys[SDLK_SPACE]) {
      camera_pos += 1 * camera_up * delta_time;
      moving = true;
    }
    if (keys[SDLK_LSHIFT]) {
      camera_pos += -1 * camera_up * delta_time;
      moving = true;
    }
    if (keys[SDLK_LEFT]) {
      cam_angle_up -= 45 * 2 * 3.14159 / 180 * delta_time;
      moving = true;
    }
    if (keys[SDLK_RIGHT]) {
      cam_angle_up += 45 * 2 * 3.14159 / 180 * delta_time;
      moving = true;
    }
    if (keys[SDLK_1]) {
      glUniform1i(u_scene, 1);
      moving = true;
    }
    if (keys[SDLK_2]) {
      glUniform1i(u_scene, 2);
      moving = true;
    }
    if (keys[SDLK_3]) {
      glUniform1i(u_scene, 3);
      moving = true;
    }
    if (keys[SDLK_4]) {
      glUniform1i(u_scene, 4);
      moving = true;
    }

    glUniform1f(u_blend, moving ? 1 : std::max(0.05, 0.1 * std::exp((1 - frames_still) / 100.0)));

    if (moving) {
      frames_still = 1;
    } else {
      frames_still += 1;
    }

    camera_forward = vec3{cos(cam_angle_up), 0, sin(cam_angle_up)};
    camera_right = cross(camera_forward, {0, 1, 0});
    camera_up = cross(camera_right, camera_forward);
    glUniform1f(g_seed, float(rand()));
    glUniform3fv(u_camera_pos, 1, camera_pos.e);
    glUniform3fv(u_camera_forward, 1, camera_forward.e);
    glUniform3fv(u_camera_up, 1, camera_up.e);
    glUniform3fv(u_camera_right, 1, camera_right.e);

    // First check to see if we should quit (from SDL)
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      switch (event.type) {
        case SDL_KEYDOWN: {
          keys[event.key.keysym.sym] = true;
        }
          break;
        case SDL_KEYUP: {
          keys[event.key.keysym.sym] = false;
        }
          break;
        case SDL_QUIT: {
          SDL_Quit();
          return 0;
        }
        break;
      }
    }
    // Set our black background
    //glClearColor(0.0, 0.0, 0.0, 1.0);
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // Now we draw the triangles. There are 3 points to draw
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    GL_DEBUG(glDrawArrays(GL_TRIANGLES, 0, 3));
    // Swap the output
    SDL_GL_SwapWindow(window);
  }
  SDL_Quit();
  return 0;
}
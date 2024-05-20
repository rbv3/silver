// Created by Marcos Oliveira <mhco@cin.ufpe.br> on 05/20/2024.
// Copyright (c)

#include <format>
#include <iostream>
#include <string_view>

#define GLAD_GL_IMPLEMENTATION
#include <glad/glad.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <glm/vec3.hpp>

#include "silver/window.h"

namespace silver {
Window::Window(int width, int height, const std::string_view& title)
    : window_(glfwCreateWindow(width, height, title.data(), nullptr, nullptr)) {
  glfwMakeContextCurrent(window_);
}

Window::~Window() {
  glfwDestroyWindow(window_);
}

void Window::MainLoop() {
  while (!glfwWindowShouldClose(window_)) {
    glfwPollEvents();

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    glClear(GL_COLOR_BUFFER_BIT);

    glfwSwapBuffers(window_);
  }
}
}  // namespace silver

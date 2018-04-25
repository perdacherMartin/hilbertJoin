
#ifndef ENERGY_CONSUMPTION_H
#define ENERGY_CONSUMPTION_H

#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <termios.h>
#include <unistd.h>
#include <sys/ioctl.h>

#include <stdio.h>
#include <stdlib.h>

#define SLEEPVAL 10000

class Hioki{
public:
    Hioki();
    ~Hioki();
    void start();
    void stop();
    void reset();
    double getWH();
private:
    const char *device = "/dev/ttyS0";
    int fd;

    int getFd();
    void setFd(int filedescriptor);
    void rs_232_config();
    void rs_232_open();

};

#endif

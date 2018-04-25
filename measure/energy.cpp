
#include "energy.h"

#include <sys/ioctl.h>

Hioki::Hioki(){
    this->rs_232_open();
    this->rs_232_config();
}

Hioki::~Hioki(){
    close(fd);
}

void Hioki::start(){
    const char *query = ":INTEGrate:STATe START\n";

    int nbytes = write(getFd(), query, strlen(query));
    if (nbytes < 0){
        perror("write of 'start' failed!\n");
        exit(1);
    }
    // printf("start: %d written\n",nbytes);
    usleep( SLEEPVAL * (strlen(query) ) );
}

int Hioki::getFd(){
    return fd;
}

void Hioki::setFd(int filedescriptor){
    fd = filedescriptor;
}

void Hioki::reset(){
    const char *query = ":INTEGrate:STATe RESET\n";
    int nbytes = write(getFd(), query, strlen(query));
    if (nbytes < 0){
        perror("write of 'reset' failed!\n");
        exit(1);
    }

    usleep(SLEEPVAL * (strlen(query) ) );
}

void Hioki::stop(){
    const char *query = ":INTEGrate:STATe STOP\n";

    int nbytes = write(getFd(), query, strlen(query));
    if (nbytes < 0){
        // printf("fd: %d\n", getFd());
        perror("write of 'stop' failed!\n");
        exit(1);
    }

    usleep(SLEEPVAL * (strlen(query) ) );

}

double Hioki::getWH(){
    char buffer[1024];
    char *bufptr;
    int nbytes,nread;
    double valueWH = 0.0;
    const char *query = ":MEASure? WH\n";

    nbytes = write(getFd(), query, strlen(query));
    if (nbytes < 0){
        perror("query of 'measure' failed!\n");
        exit(1);
    }

    usleep(SLEEPVAL * (strlen(query) + 50) );
    bufptr = buffer;
    // are there some bytes available on input
    ioctl(getFd(), FIONREAD, &nbytes);
    if ( nbytes > 0 ){
        nread = read(getFd(), bufptr, nbytes);
        bufptr[nread] = '\0';

        // receiving data in format:
        // WH +0.00053E+3
        bufptr += 3;

        valueWH = atof(bufptr);
    }else{
        fprintf(stderr, "Error in parsing value!\n");
        exit(1);
    }

    return valueWH;
}

void Hioki::rs_232_open(){
    int fdesc = open(device, O_CREAT | O_RDWR | O_NOCTTY | O_NDELAY);
    setFd(fdesc);
    // printf("fd: %d\n", getFd() );
    if(fd == -1) {
        perror("failed to open port\n" );
        exit(1);
    }else
        fcntl(getFd(), F_SETFL, 0);
}

void Hioki::rs_232_config(){
    struct termios options;
    // configuring port
    tcgetattr(getFd(), &options);
    cfsetispeed(&options, B9600); // BAUDRATE to 9600
    cfsetospeed(&options, B9600);

    options.c_cflag = (options.c_cflag & ~CSIZE) | CS8; // 8-bit chars
    options.c_iflag &= ~IGNBRK;
    options.c_lflag = 0;
    options.c_oflag = 0;                // no remapping, no delays
    options.c_cc[VMIN]  = 0;            // read doesn't block
    options.c_cc[VTIME] = 5;            // 0.5 seconds read timeout

    options.c_iflag &= ~(IXON | IXOFF | IXANY); // shut off xon/xoff ctrl

    options.c_cflag |= (CLOCAL | CREAD);// ignore modem controls,
                                    // enable reading
    options.c_cflag &= ~(PARENB | PARODD);      // shut off parity
    options.c_cflag |= 0; // 0 no parity, PARENB|PARODD odd parity, PARENB (enable parity and use even), PARENB|PARODD|CMSPAR (mark parity), and PARENB|CMSPAR (space parity).
    options.c_cflag &= ~CSTOPB;
    options.c_cflag &= ~CRTSCTS;

    /*
     * Enable the receiver and set local mode...
     */
    options.c_cflag |= (CLOCAL | CREAD);
    /*
     * Set the new options for the port...
     */
    if ( tcsetattr(getFd(), TCSANOW, &options) < 0 ){
        perror("Failed to apply settings\n");
        exit(1);
    }

}

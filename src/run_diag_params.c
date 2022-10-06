#include <stdlib.h>
#include <stdio.h>
#include <unistd.h> /* for fork */
#include <sys/types.h> /* for pid_t */
#include <sys/wait.h> /* for wait */
#include <fcntl.h>
#include <time.h>
#include "constants.h"

int create_process(char * cmdName, char ** argv, char * logfile);

int main(int argc, char ** argv)
{
    // char * logfile = "logs/run_diag_params.log";
    char * logfile = "/dev/null";

    char size[] = "10";
    char coupling[] = "0.0";
    char hopup[] = "1.0";
    char hopdn[] = "1.0";
    char hopping[] = "1.0";
    char numRuns[] = "100";
    int nospin = 0;
    DTYPE disorder_low = 5.0;
    DTYPE disorder_high = 20.0;
    int disorder_samples = 11;

    DTYPE disorder;
    char disorder_str[32];

    time_t start_time = time(NULL);
    struct tm dt = *localtime(&start_time);
    FILE * ofile = fopen(logfile, "a");
    if(ofile != NULL)
    {
        fprintf(ofile, "\n\nSTARTING RUN AT %d-%02d-%02d %02d:%02d:%02d\n",
                dt.tm_year + 1900, dt.tm_mon + 1, dt.tm_mday, dt.tm_hour,
                dt.tm_min, dt.tm_sec);
        fclose(ofile);
    }

    int i;
    DTYPE coeff = (disorder_high - disorder_low) / (DTYPE) (disorder_samples-1);
    for(i = 0; i < disorder_samples; i++)
    {
        disorder = disorder_low + coeff * i;
        sprintf(disorder_str, "%lf", disorder);

        printf("DIAG: W = %6.3f... ", disorder);
        fflush(stdout);
        start_time = time(NULL);
        
        if(nospin == 0)
        {
            char * cmdName = "build/exact_diag_simulation";
            char * cmdArgv[] = {"exact_diag_simulation",
                                "-s", size,
                                "-c", coupling,
                                "-u", hopup,
                                "-d", hopdn,
                                "-n", numRuns,
                                "-w", disorder_str,
                                NULL};
            create_process(cmdName, cmdArgv, logfile);
        }
        else
        {
            char * cmdName = "build/exact_diag_simulation";
            char * cmdArgv[] = {"exact_diag_simulation",
                                "-s", size,
                                "-t", hopping,
                                "-n", numRuns,
                                "-w", disorder_str,
                                "-p",
                                NULL};
            create_process(cmdName, cmdArgv, logfile);
        }

        printf("Done in %lds\n", (time(NULL) - start_time));
    }
    // wait();

    return(0);
}

int create_process(char * cmdName, char ** argv, char * logfile)
{
    // Spawn a child process
    pid_t pid = fork();

    // Forking creates two copies
    // of this very program.
    // Is the copy that is running
    // the parent or the child?

    if(pid == 0) // This is the child process
    {
        /* open file for writing */
        int fd = open(logfile, O_WRONLY | O_APPEND | O_CREAT, S_IRUSR | S_IWUSR);

        dup2(fd, 1);    /* make stdout a copy of fd (> /dev/null) */
        // dup2(fd, 2);    /* ...and same with stderr */
        close(fd);      /* close fd */

        /* stdout and stderr now write to /dev/null */

        execv(cmdName, argv);
        // if execv() fails then
        // the next statement will be
        // executed.
        exit(127);
    }
    else // This is the parent process
    {
        // Wait for the child to finish
        waitpid(pid, 0, 0);
    }
    return(0);
}


#pragma once
#include <string>
#include <string.h>

class option_handler {
 public:
    option_handler(int argc, char **argv) {
        m_argc = argc;
        m_argv = argv;
    }

    bool get_value(const char* option, int& value) {
        int index;
        bool match = get_match(option, index);
        if(match) {
            value = atoi(m_argv[index + 1]);
        }
        return match;
    }

    bool get_value(const char* option, double& value) {
        int index;
        bool match = get_match(option, index);
        if(match) {
            value = std::stod(m_argv[index + 1]);
        }
        return match;
    }

    bool get_value(const char* option, bool& value) {
        int value_int;
        bool match = get_value(option, value_int);
        if(match) {
            value = (bool) value_int;
        }
        return match;
    }
    
    bool get_value(const char* option, char *value){
    	  int index;
        bool match = get_match(option, index);
        if(match) {
            sprintf(value, "%s", m_argv[index + 1]);
        }
        return match;
    }

 private:

    bool get_match(const char* option, int& index)
    {
        for(index = 0; index < m_argc; index++) {
            char *arg_i = m_argv[index];
            if(!strcmp(option, arg_i)) return true;
        }
        return false;
    }

    int m_argc;
    char **m_argv;
};

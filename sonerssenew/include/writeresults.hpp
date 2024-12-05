#ifndef _WRITERESULTS_HPP_DEFINED_
#define _WRITERESULTS_HPP_DEFINED_

class writeresults {
public:
    //template <typename... Args>
    //void save_bin_data(const char* , const Args&... );

    void save_bin_data(const char* , const char** headers, const double* values, int no_of_observables);
    
    //private:
    //    template <typename T>
    //    void writeHeader(std::ostream& file, const T& arg);

    //    template <typename T, typename... Args>
    //    void writeHeader(std::ostream& file, const T& arg, const Args&... args);
};

#endif // WRITERESULTS_H

#pragma once

template<class T> constexpr T get_numeric_limit(){ return 1e-13;
    /*
                                            cout<<"ef: "<<std::numeric_limits<float>::max_digits10<<endl;
    cout<<"ed: "<<std::numeric_limits<double>::max_digits10<<endl;
    cout<<"el: "<<std::numeric_limits<long double>::max_digits10<<endl;
*/
}
template<> constexpr float get_numeric_limit<float>(){
    // abs limit is 9 digits
    return 1e-7;
}
template<> constexpr double get_numeric_limit<double>(){
    // abs limit is 17 digits
    return 1e-13;
}
template<> constexpr long double get_numeric_limit<long double>(){
    // abs limit is 21 digits
    return 1e-15 ;
}

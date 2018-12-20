
template<class Source, class Messanger>
auto sendSource(Source src, Messanger&& message) -> void {
    while(src.has()) {
        message(src.get()); 
        src.advance(); 
    }
}

template<class Iter>
struct IteratedSource {
    Iter beginning;
    Iter ending; 
    auto has() const noexcept(noexcept(beginning != ending)) -> bool {
        return beginning != ending; 
    }
    auto get() noexcept(noexcept(*beginning)) -> decltype(*beginning) {
        return *beginning; 
    }
    auto get() const noexcept(noexcept(*beginning)) -> decltype(*beginning) {
        return *beginning; 
    }
    auto advance() noexcept(noexcept(++beginning)) -> void {
        ++beginning; 
    }
};

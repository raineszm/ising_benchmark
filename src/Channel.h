#pragma once

#include <queue>
#include <mutex>
#include <condition_variable>

template <class T>
class Channel {
    private:
        std::mutex mtx;
        std::condition_variable cv;
        std::queue<T> queue;

        T _take() {
            T item = queue.front();
            queue.pop();
            return item;
        }

    public:

        void put(const T& item) {
            {
                std::lock_guard<std::mutex> lock(mtx);
                queue.push(item);
            }
            cv.notify_one();
        }

        T take() {
            std::unique_lock<std::mutex> lock(mtx);

            if (queue.empty())
                cv.wait(lock, [this]() { return !queue.empty(); });

            return _take();
        }
};

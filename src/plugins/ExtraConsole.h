//
// Created by alex on 4/4/25.
//

#pragma once

// Logger
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <cstdlib>
#include <cstdio>

// Plotter
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <type_traits>

namespace Terminal {
    class Logger {
    public:
        explicit Logger(const std::string& fifo_path = "/tmp/plugin_fifo")
                : fifo_path_(fifo_path + "_" + std::to_string(next_id_++)), fifo_stream_() {
            setup_fifo();
            launch_terminal();
            open_stream();
        }

        ~Logger() {
            if (fifo_stream_.is_open()) {
                fifo_stream_.close();
            }
            unlink(fifo_path_.c_str());
        }

        std::ostream& stream() {
            return fifo_stream_;
        }

    private:
        std::string fifo_path_;
        std::ofstream fifo_stream_;
        inline static int next_id_ = 0;

        void setup_fifo() {
            // Create a named FIFO if it doesn't exist
            if (access(fifo_path_.c_str(), F_OK) == -1) {
                if (mkfifo(fifo_path_.c_str(), 0666) != 0) {
                    perror("mkfifo failed");
                    std::exit(EXIT_FAILURE);
                }
            }
        }

        void launch_terminal() {
            // Launch gnome-terminal and run cat on the FIFO
            std::string cmd = "gnome-terminal -- bash -c \"cat " + fifo_path_ + "; exec bash\"";
            int ret = std::system(cmd.c_str());
            if (ret != 0) {
                std::cerr << "Failed to launch gnome-terminal" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            // Optional: brief sleep to ensure terminal is ready
            usleep(300'000);
        }

        void open_stream() {
            fifo_stream_.open(fifo_path_, std::ios::out);
            if (!fifo_stream_.is_open()) {
                std::cerr << "Failed to open FIFO for writing\n";
                std::exit(EXIT_FAILURE);
            }
        }
    };

    template<typename T, typename check = typename std::enable_if_t<std::is_arithmetic_v<T>>>
    class Plotter {
    public:
        explicit Plotter(const std::string& name = "", size_t width = 80, size_t height = 15)
        : name_(name), width_(width), height_(height), data_(width, 0) {}

        void add_value(const T& value) {
            data_.push_back(value);
            if (data_.size() > width_) {
                data_.erase(data_.begin());
            }
        }

        void set_values(const std::vector<T>& values) {
            // do sub sampling
            if (values.size() > width_) {
                const double step = static_cast<double>(values.size()) / static_cast<double>(width_);
                for (int idx = 0; idx < width_; idx++) {
                    data_[idx] = values[static_cast<int>(static_cast<double>(idx) * step)];
                }
            }
            // just copy over
            else std::copy(values.begin(), values.end(), data_.begin());
        }

        [[nodiscard]] std::string get_plot_string() const {
            // Determine min and max
            T min_val = *std::min_element(data_.begin(), data_.end());
            T max_val = *std::max_element(data_.begin(), data_.end());
            return get_plot_string(min_val, max_val);
        }

        [[nodiscard]] std::string get_plot_string(const T& min_val, const T& max_val) const {
            T range = std::max(max_val - min_val, 1e-5); // prevent div-by-zero

            std::ostringstream out;
            out << "\033[H\033[J"; // ANSI: home and clear
            out << name_ << "\n";

            // Top-down row rendering
            for (int row = height_ - 1; row >= 0; --row) {
                T threshold = min_val + range * row / height_;
                for (T val : data_) {
                    T clamped = std::clamp(val, min_val, max_val);
                    out << (clamped >= threshold ? "█" : " ");
                }
                out << "\n";
            }

            out << std::string(width_, '-') << "\n";
            out << std::fixed << std::setprecision(2);
            out << "Y Range: [" << min_val << ", " << max_val << "]\n";
            return out.str();
        }

    private:
        std::string name_;
        size_t width_, height_;
        std::vector<T> data_;
    };

    template<typename T>
    class PlotLogger {
    public:
        PlotLogger(const std::string& name) : logger_("/tmp/ls1mardyn_plugin_fifo"), plotter_(name) {}

        void add(const T& value) {
            plotter_.add_value(value);
            logger_.stream() << plotter_.get_plot_string() << std::flush;
        }

        void set(const std::vector<T>& values) {
            plotter_.set_values(values);
            logger_.stream() << plotter_.get_plot_string() << std::flush;
        }

        void set(const std::vector<T>& values, const T& min, const T& max) {
            plotter_.set_values(values);
            logger_.stream() << plotter_.get_plot_string(min, max) << std::flush;
        }

    private:
        Logger logger_;
        Plotter<T> plotter_;
    };
}


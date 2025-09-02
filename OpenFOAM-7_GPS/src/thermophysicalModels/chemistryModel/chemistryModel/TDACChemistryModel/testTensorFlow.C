#include <iostream>
#include <tensorflow/core/framework/tensor.h>
#include <tensorflow/core/public/session.h>
#include <tensorflow/core/platform/env.h>

int main() {
    // Initialize the TensorFlow session
    tensorflow::Session* session;
    tensorflow::SessionOptions session_options;
    tensorflow::Status status = tensorflow::NewSession(session_options, &session);
    if (!status.ok()) {
        std::cerr << "Error creating TensorFlow session: " << status.ToString() << std::endl;
        return -1;
    }

    // Load the SavedModel
    tensorflow::GraphDef graph_def;
    status = tensorflow::ReadBinaryProto(tensorflow::Env::Default(), "/path/to/saved_model/saved_model.pb", &graph_def);
    if (!status.ok()) {
        std::cerr << "Error reading SavedModel: " << status.ToString() << std::endl;
        return -1;
    }

    status = session->Create(graph_def);
    if (!status.ok()) {
        std::cerr << "Error creating session with SavedModel: " << status.ToString() << std::endl;
        return -1;
    }

    // Prepare input data (example)
    std::vector<float> input_data = {1.0, 2.0, 3.0}; // Example input data

    // Create input tensor
    tensorflow::Tensor input_tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, input_data.size()}));
    auto input_tensor_mapped = input_tensor.tensor<float, 2>();
    for (int i = 0; i < input_data.size(); ++i) {
        input_tensor_mapped(0, i) = input_data[i];
    }

    // Prepare output tensor
    std::vector<tensorflow::Tensor> output_tensors;

    // Run the session to perform inference
    status = session->Run({{ "input_tensor_name", input_tensor }}, { "output_tensor_name" }, {}, &output_tensors);
    if (!status.ok()) {
        std::cerr << "Error running session: " << status.ToString() << std::endl;
        return -1;
    }

    // Get the output tensor value
    tensorflow::Tensor output_tensor = output_tensors[0];
    auto output_tensor_mapped = output_tensor.tensor<float, 2>();

    // Print the prediction result
    std::cout << "Prediction result:" << std::endl;
    for (int i = 0; i < output_tensor_mapped.dimension(1); ++i) {
        std::cout << output_tensor_mapped(0, i) << std::endl;
    }

    // Cleanup
    status = session->Close();
    if (!status.ok()) {
        std::cerr << "Error closing TensorFlow session: " << status.ToString() << std::endl;
        return -1;
    }

    return 0;
}

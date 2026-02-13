package org.gelfgren;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

/**
 * Main entry point for Gelfgren library.
 * Handles native library loading.
 */
public class Gelfgren {
    private static boolean loaded = false;

    static {
        loadNativeLibrary();
    }

    /**
     * Load the native Gelfgren library.
     * Attempts to load from java.library.path first, then extracts from JAR.
     */
    private static synchronized void loadNativeLibrary() {
        if (loaded) {
            return;
        }

        // Try loading from java.library.path
        try {
            System.loadLibrary("gelfgren");
            loaded = true;
            return;
        } catch (UnsatisfiedLinkError e) {
            // Will try extracting from JAR
        }

        // Determine library name based on OS
        String osName = System.getProperty("os.name").toLowerCase();
        String libName;
        if (osName.contains("win")) {
            libName = "gelfgren.dll";
        } else if (osName.contains("mac")) {
            libName = "libgelfgren.dylib";
        } else {
            libName = "libgelfgren.so";
        }

        // Try to extract and load from resources
        try {
            InputStream is = Gelfgren.class.getResourceAsStream("/native/" + libName);
            if (is == null) {
                throw new UnsatisfiedLinkError(
                    "Native library not found in JAR: /native/" + libName);
            }

            Path tempFile = Files.createTempFile("gelfgren", libName);
            tempFile.toFile().deleteOnExit();
            Files.copy(is, tempFile, StandardCopyOption.REPLACE_EXISTING);
            is.close();

            System.load(tempFile.toAbsolutePath().toString());
            loaded = true;
        } catch (IOException e) {
            throw new UnsatisfiedLinkError(
                "Failed to load native library: " + e.getMessage());
        }
    }

    /**
     * Get the library version.
     */
    public static String getVersion() {
        return "0.1.0";
    }

    /**
     * Get the last error message from the native library.
     */
    public static native String getLastErrorMessage();

    /**
     * Clear the last error message.
     */
    public static native void clearLastError();
}

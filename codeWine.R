library(shiny)
library(ggplot2)

ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            sliderInput("seeds", "Seeded Grapes", 0, 10000, 5000, step = 100),
            sliderInput("survival", "Survival Rate", 0, 1, 0.8, step = 0.01),
            sliderInput("harvest", "Harvest Rate", 0, 1, 0.5, step = 0.01),
            sliderInput("quality", "Quality Rate", 0, 1, 0.6, step = 0.01),
            plotOutput("compartment_plot")
        ),
        mainPanel(
            tableOutput("compartment_table")
        )
    )
)

server <- function(input, output) {
    
    # Define compartmental equations
    compartment_equations <- function(time, state, parameters) {
        
        # Assign initial values to compartments
        seeds <- state[1]
        survived <- state[2]
        harvested <- state[3]
        wine <- state[4]
        discarded <- state[5]
        
        # Extract parameter values
        sowing_rate <- parameters[1]
        survival_rate <- parameters[2]
        harvest_rate <- parameters[3]
        winery_yield <- parameters[4]
        
        # Calculate the rates of change
        dsdt <- -sowing_rate * seeds
        dsurvived_dt <- sowing_rate * seeds - survival_rate * survived
        dharvested_dt <- survival_rate * survived - harvest_rate * harvested
        dwine_dt <- harvest_rate * harvested * winery_yield
        ddiscarded_dt <- sowing_rate * seeds * (1 - survival_rate)
        
        # Return rates of change as a list
        list(c(dsdt, dsurvived_dt, dharvested_dt, dwine_dt, ddiscarded_dt))
    }
    
    # Define initial state and time vector
    initial_state <- c(seeds = 1000, survived = 0, harvested = 0, wine = 0, discarded = 0)
    times <- seq(0, 365, by = 1)
    
    # Define parameter vector and variable names
    initial_parameters <- c(sowing_rate = 0.002, survival_rate = 0.8, harvest_rate = 0.01, winery_yield = 0.8)
    parameter_names <- c("Sowing Rate", "Survival Rate", "Harvest Rate", "Winery Yield")
    
    # Create reactive values for parameters and state
    parameters <- reactiveValues(values = initial_parameters)
    state <- reactiveValues(values = initial_state)
    
    # Update state based on parameter values
    observeEvent(parameters$values, {
        ode_out <- ode(y = state$values, times = times, func = compartment_equations, parms = parameters$values)
        state$values <- ode_out[ ,5]
    })
    
    # Update plot based on state
    output$plot <- renderPlot({
        plot(times, state$values[,2], type = "l", xlab = "Time (days)", ylab = "Population", main = "Grapevine Population")
        lines(times, state$values[,3], col = "red")
        lines(times, state$values[,4], col = "blue")
        legend("topright", legend = c("Survived", "Harvested", "Wine"), col = c("black", "red", "blue"), lty = 1)
    })
    
    # Update parameter values based on sliders
    observe({
        parameters$values[1] <- input$sowing_rate / 1000
        parameters$values[2] <- input$survival_rate / 100
        parameters$values[3] <- input$harvest_rate / 1000
        parameters$values[4] <- input$winery_yield / 100
    })
    
    # Compartmental equations
    ode_fun <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
            s <- seeds
            g <- germinated
            v <- vines
            f <- fruit
            h <- harvested
            w <- wine
            d <- discarded
            
            ds_dt <- -s * v * r_sow
            dg_dt <- s * v * r_sow - g * r_survival
            dv_dt <- g * r_survival - f * r_vine_to_fruit - h * r_vine_to_harvest - d * r_vine_to_discard
            df_dt <- f * r_vine_to_fruit - h * r_fruit_to_harvest - d * r_fruit_to_discard
            dh_dt <- h * r_harvest_to_wine * r_winery_yield - w * r_wine_to_bottle
            dw_dt <- h * r_harvest_to_wine * (1 - r_winery_yield) + w * r_wine_to_bottle
            dd_dt <- d * r_vine_to_discard + f * r_fruit_to_discard
            return(list(c(ds_dt, dg_dt, dv_dt, df_dt, dh_dt, dw_dt, dd_dt)))
        })
    }
    
    # Run the simulation
    output$plot <- renderPlot({
        time <- seq(0, input$duration, by = 0.1)
        state <- as.numeric(c(seeds = input$seeds, germinated = 0, vines = 0, fruit = 0, harvested = 0, wine = 0, discarded = 0))
        parameters <- reactiveValuesToList(parameters)
        rates <- reactive({
            ode(y = state, times = time, func = ode_fun, parms = parameters)
        })
        rates_df <- as.data.frame(rates())
        rates_df$time <- time
        rates_df_long <- gather(rates_df, key = "compartment", value = "count", -time)
        ggplot(data = rates_df_long, aes(x = time, y = count, group = compartment, color = compartment)) +
            geom_line() +
            xlab("Time (months)") +
            ylab("Number of grapes") +
            scale_color_discrete(name = "Compartment")
    })
    
    # Sliders
    observe({
        sliderInput("sowing_rate", "Sowing rate (kg/ha)", min = 10, max = 500, value = 100)
        sliderInput("survival_rate", "Survival rate (%)", min = 0, max = 100, value = 50)
        sliderInput("harvest_rate", "Harvest rate (kg/ha)", min = 100, max = 5000, value = 1000)
        sliderInput("winery_yield", "Winery yield (%)", min = 10, max = 90, value = 70)
    })
    
    # Parameter table
    output$params_table <- renderTable({
        parameters <- parameters$values
        data.frame(Parameter = c("Sowing rate (kg/ha)", "Survival rate (%)", "Harvest rate (kg/ha)", "Winery yield (%)"), Value = parameters,
                   stringsAsFactors = FALSE)
    })
    
    
    output$plot <- renderPlot({
        # Solve differential equations
        out <- ode(y = init, times = times, func = grape_model, parms = parameters)
        
        # Plot compartments over time
        df <- data.frame(Time = out[,1], Seeded = out[,2], Survived = out[,3], Harvested = out[,4], Wine = out[,5], Discarded = out[,6])
        df_long <- gather(df, Compartment, Population, -Time)
        ggplot(df_long, aes(x = Time, y = Population, color = Compartment)) +
            geom_line() +
            theme_minimal() +
            labs(x = "Time (months)", y = "Population (kg/ha)", color = "") +
            scale_color_manual(values = c("#8DA0CB", "#FC8D62", "#E78AC3", "#A6D854", "#66C2A5"))
    })
    
    
    output$slider_sowing_rate <- renderUI({
        sliderInput("sowing_rate", "Sowing rate (kg/ha)", min = 0, max = 1000, value = 250, step = 1)
    })
    
    output$slider_survival_rate <- renderUI({
        sliderInput("survival_rate", "Survival rate (%)", min = 0, max = 100, value = 50, step = 1)
    })
    
    output$slider_harvest_rate <- renderUI({
        sliderInput("harvest_rate", "Harvest rate (kg/ha)", min = 0, max = 1000, value = 250, step = 1)
    })
    
    output$slider_winery_yield <- renderUI({
        sliderInput("winery_yield", "Winery yield (%)", min = 0, max = 100, value = 75, step = 1)
    })
    
}

shinyApp(ui = ui, server = server)
